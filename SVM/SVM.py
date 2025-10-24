import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter
import joblib
import json
import os

from sklearn.model_selection import GridSearchCV, learning_curve, validation_curve
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import (classification_report, confusion_matrix, roc_curve, auc,
                             precision_recall_curve, average_precision_score, matthews_corrcoef)

# --- Global Constants ---
# Amino acid and physicochemical properties
AMINO_ACIDS = sorted(list("ACDEFGHIKLMNPQRSTVWY"))
N_TERMINUS_LEN = 22
HYDROPHOBICITY_SCALE_KD = {
    'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5, 'Q': -3.5, 'E': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5, 'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6, 'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2
}
TM_HELIX_PROPENSITY_SCALE = {
    'A': 0.25, 'R': -1.80, 'N': -0.64, 'D': -0.72, 'C': 0.04, 'Q': -0.69, 'E': -0.62, 'G': 0.16, 'H': -0.40, 'I': 0.73, 'L': 0.53, 'K': -1.10, 'M': 0.26, 'F': 0.61, 'P': -0.07, 'S': -0.26, 'T': -0.18, 'W': 0.37, 'Y': 0.02, 'V': 0.54
}
ALPHA_HELIX_PROPENSITY = {
    'A': 1.42, 'R': 0.98, 'N': 0.67, 'D': 1.01, 'C': 0.70, 'Q': 1.11, 'E': 1.51, 'G': 0.57, 'H': 1.00, 'I': 1.08, 'L': 1.21, 'K': 1.16, 'M': 1.45, 'F': 1.13, 'P': 0.57, 'S': 0.77, 'T': 0.83, 'W': 1.08, 'Y': 0.69, 'V': 1.06
}
CHARGE_SCALE = {'D': -1, 'E': -1, 'K': 1, 'R': 1}

# --- Feature Engineering Functions ---

def calculate_average_property(sequence, scale):
    """Calculates the average value of a physicochemical property for a sequence."""
    if not sequence: return 0.0
    return sum(scale.get(aa, 0.0) for aa in sequence) / len(sequence)

def create_augmented_features(sequence):
    """Generates a feature vector from a protein sequence."""
    if not isinstance(sequence, str): sequence = ""
    n_term_seq = sequence[:N_TERMINUS_LEN].upper()
    n_region, h_region = n_term_seq[0:5], n_term_seq[5:16]
    
    if not n_term_seq: 
        # Return a zero vector if the sequence is empty
        return [0.0] * (len(AMINO_ACIDS) + 6)
        
    counts, total_len = Counter(n_term_seq), len(n_term_seq)
    
    aa_composition = [(counts.get(aa, 0) / total_len) for aa in AMINO_ACIDS]
    n_region_charge = sum(CHARGE_SCALE.get(aa, 0) for aa in n_region)
    h_region_hydrophobicity_kd = calculate_average_property(h_region, HYDROPHOBICITY_SCALE_KD)
    h_region_propensity_tm = calculate_average_property(h_region, TM_HELIX_PROPENSITY_SCALE)
    avg_hydrophobicity_total = calculate_average_property(n_term_seq, HYDROPHOBICITY_SCALE_KD)
    avg_propensity_tm_total = calculate_average_property(n_term_seq, TM_HELIX_PROPENSITY_SCALE)
    avg_propensity_alpha_total = calculate_average_property(n_term_seq, ALPHA_HELIX_PROPENSITY)
    
    return aa_composition + [
        n_region_charge, h_region_hydrophobicity_kd, h_region_propensity_tm,
        avg_hydrophobicity_total, avg_propensity_tm_total, avg_propensity_alpha_total
    ]

# --- Main Pipeline Function ---

def main():
    """Main function to run the complete ML pipeline."""
    
    # Set plot style for consistency
    sns.set_theme(context="notebook", style="darkgrid", palette="muted")

    # --- Create a directory for plots ---
    plots_dir = "plots"
    if not os.path.exists(plots_dir):
        os.makedirs(plots_dir)
        print(f"Created directory: {plots_dir}")

    # --- 1. Data Loading and Preparation ---
    print("--- 1. Loading and Preparing Data ---")
    try:
        df = pd.read_csv("merged_dataset_with_seqs.tsv", sep="\t")
        print("Dataset loaded successfully.")
    except FileNotFoundError:
        print("Error: 'merged_dataset_with_seqs.tsv' not found. Please place it in the same directory.")
        exit()

    train_df = df[df["Type"] == "train"].copy()
    test_df = df[df["Type"] == "test"].copy()
    print(f"Training set size: {len(train_df)}")
    print(f"Test set size: {len(test_df)}")

    # --- 2. Feature Engineering (Augmented Features) ---
    print("\n--- 2. Engineering Augmented Features ---")
    
    X_train_list = train_df['Sequence'].apply(create_augmented_features).tolist()
    y_train = train_df['label'].values
    X_train = np.array(X_train_list)
    
    feature_names = ([f"AA_{aa}_freq" for aa in AMINO_ACIDS] + 
                     ["n_region_charge", "h_region_hydrophobicity", "h_region_TM_propensity", 
                      "total_avg_hydrophobicity", "total_avg_TM_propensity", "total_avg_alpha_propensity"])
    
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    print(f"Training feature engineering complete. Shape: {X_train_scaled.shape}")
    
    if not test_df.empty:
        X_test_list = test_df['Sequence'].apply(create_augmented_features).tolist()
        y_test = test_df['label'].values
        X_test = np.array(X_test_list)
        X_test_scaled = scaler.transform(X_test)
        print(f"Test feature engineering complete. Shape: {X_test_scaled.shape}")
    else:
        X_test_scaled, y_test = np.array([]), np.array([])

    # --- 3. Feature Selection using RandomForest ---
    print("\n--- 3. Running Feature Selection with RandomForest ---")
    rf_selector = RandomForestClassifier(n_estimators=150, random_state=42, n_jobs=-1)
    rf_selector.fit(X_train_scaled, y_train)
    importances = rf_selector.feature_importances_
    feature_importance_df = pd.DataFrame({'Feature': feature_names, 'Importance': importances}).sort_values(by='Importance', ascending=False)
    
    plt.figure(figsize=(12, 10))
    sns.barplot(x='Importance', y='Feature', data=feature_importance_df, hue='Feature', palette='viridis_r', legend=False)
    plt.title('Augmented Feature Importances from RandomForest', fontsize=16)
    plt.tight_layout()
    plt.savefig(os.path.join(plots_dir, "svm_augmented_feature_importances.png"), dpi=300)
    # plt.show() # Disabled for non-interactive script execution
    
    print("Top 15 most important features:\n", feature_importance_df.head(15).to_string())

    # --- 3.5. Feature Filtering ---
    print("\n--- 3.5. Filtering Dataset to Top 15 Features ---")
    top_15_feature_names = feature_importance_df.head(15)['Feature'].tolist()
    top_15_indices = [feature_names.index(feature) for feature in top_15_feature_names]
    
    X_train_selected = X_train_scaled[:, top_15_indices]
    if X_test_scaled.shape[0] > 0:
        X_test_selected = X_test_scaled[:, top_15_indices]
    else:
        X_test_selected = np.array([])
        
    print(f"Filtered training features shape: {X_train_selected.shape}")
    print(f"Filtered test features shape: {X_test_selected.shape}")
    
    # --- 4. SVM Hyperparameter Optimization with GridSearchCV (with Caching) ---
    print("\n--- 4. Optimizing SVM Hyperparameters with GridSearchCV ---")
    best_params_filename = 'best_svm_params.json'
    
    if os.path.exists(best_params_filename):
        print(f"Found existing parameters in '{best_params_filename}'. Loading them.")
        with open(best_params_filename, 'r') as f:
            best_params = json.load(f)
        
        print(f"Loaded parameters: {best_params}")
        best_svm = SVC(**best_params, probability=True, random_state=42)
        best_svm.fit(X_train_selected, y_train)
        print("Model training complete.")
    else:
        print(f"'{best_params_filename}' not found. Running GridSearchCV.")
        param_grid = [
            {'kernel': ['linear'], 'C': [0.1, 1, 10, 100]},
            {'kernel': ['rbf'], 'C': [0.1, 1, 10, 100], 'gamma': ['scale', 'auto', 0.01, 0.1]},
            {'kernel': ['poly'], 'C': [0.1, 1, 10, 100], 'degree': [2, 3], 'gamma': ['scale', 'auto']}
        ]
        
        cv_iterator = []
        index_mapper = {orig_idx: arr_idx for arr_idx, orig_idx in enumerate(train_df.index)}
        for i in range(5):
            train_fold_indices = train_df[train_df['fold'] != i].index
            val_fold_indices = train_df[train_df['fold'] == i].index
            train_indices = [index_mapper[idx] for idx in train_fold_indices]
            val_indices = [index_mapper[idx] for idx in val_fold_indices]
            cv_iterator.append((train_indices, val_indices))
            
        grid_search = GridSearchCV(estimator=SVC(probability=True, random_state=42), 
                                   param_grid=param_grid, cv=cv_iterator, 
                                   scoring='f1', n_jobs=-1, verbose=2)
        grid_search.fit(X_train_selected, y_train)
        
        best_svm = grid_search.best_estimator_
        best_params = grid_search.best_params_
        print("\nGridSearchCV complete.")
        print(f"Best parameters found: {best_params}")
        print(f"Best cross-validation F1-score: {grid_search.best_score_:.4f}")
        
        print(f"\n--- Saving Best Parameters to {best_params_filename} ---")
        with open(best_params_filename, 'w') as f:
            json.dump(best_params, f, indent=4)
        print("Successfully saved the best parameters.")

    # --- 5. Final Model Evaluation on the Test Set ---
    print("\n--- 5. Evaluating the Best Model on the Hold-out Test Set ---")
    if X_test_selected.shape[0] > 0:
        y_pred = best_svm.predict(X_test_selected)
        y_pred_proba = best_svm.predict_proba(X_test_selected)[:, 1]

        mcc = matthews_corrcoef(y_test, y_pred)
        print(f"\nMatthews Correlation Coefficient (MCC): {mcc:.4f}")
        print("\nClassification Report on Test Set:")
        print(classification_report(y_test, y_pred, target_names=['Negative (0)', 'Positive (1)']))

        cm = confusion_matrix(y_test, y_pred)
        plt.figure(figsize=(8, 6))
        sns.heatmap(cm, annot=True, fmt='d', cmap='Blues', cbar=False,
                    xticklabels=['Negative', 'Positive'], yticklabels=['Negative', 'Positive'])
        plt.title('Confusion Matrix on Test Set', fontsize=16)
        plt.xlabel('Predicted Label', fontsize=12)
        plt.ylabel('True Label', fontsize=12)
        plt.savefig(os.path.join(plots_dir, "svm_confusion_matrix.png"), dpi=300)
        # plt.show()
        
        fpr, tpr, _ = roc_curve(y_test, y_pred_proba)
        roc_auc = auc(fpr, tpr)
        plt.figure(figsize=(8, 6))
        plt.plot(fpr, tpr, color='darkorange', lw=2, label=f'ROC curve (AUC = {roc_auc:.2f})')
        plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
        plt.grid(True)
        plt.xlabel('False Positive Rate')
        plt.ylabel('True Positive Rate')
        plt.title('Receiver Operating Characteristic (ROC) Curve', fontsize=14)
        plt.legend(loc="lower right")
        plt.savefig(os.path.join(plots_dir, "svm_roc_curve.png"), dpi=300)
        # plt.show()
        
        precision, recall, _ = precision_recall_curve(y_test, y_pred_proba)
        avg_precision = average_precision_score(y_test, y_pred_proba)
        plt.figure(figsize=(8, 6))
        plt.step(recall, precision, where='post', color='b', alpha=0.7, label=f'AP = {avg_precision:.2f}')
        plt.fill_between(recall, precision, step='post', alpha=0.3, color='b')
        plt.xlabel('Recall')
        plt.ylabel('Precision')
        plt.ylim([0.0, 1.05])
        plt.xlim([0.0, 1.0])
        plt.grid(True)
        plt.title('Precision-Recall Curve', fontsize=14)
        plt.legend(loc="upper right")
        plt.savefig(os.path.join(plots_dir, "svm_precision_recall_curve.png"), dpi=300)
        # plt.show()

    else:
        print("\nTest set is empty. Skipping final evaluation.")

    # --- 6. Additional Analyses and Plots ---
    print("\n--- 6. Performing Additional Analyses ---")

    # 6.1 MCC vs. Number of Features
    print("\n--- 6.1 Analyzing MCC vs. Number of Features ---")
    mcc_scores = []
    # Use the indices of the already sorted importance dataframe
    sorted_feature_indices = [feature_names.index(f) for f in feature_importance_df['Feature']]
    
    n_features_range = list(range(1, X_train_scaled.shape[1] + 1))
    
    for n in n_features_range:
        top_n_indices = sorted_feature_indices[:n]
        X_train_subset = X_train_scaled[:, top_n_indices]
        X_test_subset = X_test_scaled[:, top_n_indices]
        
        temp_svm = SVC(**best_params, random_state=42)
        temp_svm.fit(X_train_subset, y_train)
        y_pred_subset = temp_svm.predict(X_test_subset)
        mcc_scores.append(matthews_corrcoef(y_test, y_pred_subset))
    
    plt.figure(figsize=(10, 6))
    plt.plot(n_features_range, mcc_scores, marker='.', linestyle='-')
    plt.title('Model Performance (MCC) vs. Number of Top Features', fontsize=14)
    plt.xlabel('Number of Top Features Included', fontsize=12)
    plt.ylabel('Matthews Correlation Coefficient (MCC)', fontsize=12)
    plt.grid(True)
    plt.axvline(x=15, color='r', linestyle='--', label='15 Features (Selected)')
    plt.legend()
    plt.savefig(os.path.join(plots_dir, "svm_mcc_vs_features.png"), dpi=300)
    # plt.show()

    # 6.2 Learning Curve
    print("\n--- 6.2 Generating Learning Curve ---")
    train_sizes, train_scores, val_scores = learning_curve(
        best_svm, X_train_selected, y_train, cv=5, n_jobs=-1,
        train_sizes=np.linspace(0.1, 1.0, 10), scoring='matthews_corrcoef'
    )
    train_scores_mean = np.mean(train_scores, axis=1)
    val_scores_mean = np.mean(val_scores, axis=1)
    plt.figure(figsize=(10, 6))
    plt.plot(train_sizes, train_scores_mean, 'o-', color="r", label="Training score")
    plt.plot(train_sizes, val_scores_mean, 'o-', color="g", label="Cross-validation score")
    plt.title("Learning Curve (SVM)", fontsize=14)
    plt.xlabel("Training Examples", fontsize=12)
    plt.ylabel("MCC Score", fontsize=12)
    plt.legend(loc="best")
    plt.grid(True)
    plt.savefig(os.path.join(plots_dir, "svm_learning_curve.png"), dpi=300)
    # plt.show()
    
    # 6.3 Validation Curve for SVM 'C' parameter
    print("\n--- 6.3 Generating Validation Curve for SVM 'C' parameter ---")
    if best_params.get('kernel') != 'poly':
        param_range = np.logspace(-3, 3, 7)
        train_scores, val_scores = validation_curve(
            SVC(kernel=best_params.get('kernel', 'rbf'), gamma=best_params.get('gamma', 'scale'), random_state=42),
            X_train_selected, y_train, param_name="C", param_range=param_range,
            cv=5, scoring="accuracy", n_jobs=-1
        )
        train_scores_mean = np.mean(train_scores, axis=1)
        val_scores_mean = np.mean(val_scores, axis=1)
        plt.figure(figsize=(10, 6))
        plt.semilogx(param_range, train_scores_mean, label="Training score", color="darkorange", marker='o')
        plt.semilogx(param_range, val_scores_mean, label="Cross-validation score", color="navy", marker='o')
        plt.title("Validation Curve for SVM (C parameter)", fontsize=14)
        plt.xlabel("C (Regularization Parameter)", fontsize=12)
        plt.ylabel("Accuracy Score", fontsize=12)
        plt.legend(loc="best")
        plt.grid(True)
        plt.savefig(os.path.join(plots_dir, "svm_validation_curve.png"), dpi=300)
        # plt.show()
    else:
        print("Skipping validation curve for 'poly' kernel to keep the plot simple (2D).")
        
    # --- 7. Save the Final Model and Scaler ---
    print("\n--- 7. Saving the final trained model and scaler to disk ---")
    model_filename = 'svm_signal_peptide_model.joblib'
    scaler_filename = 'feature_scaler.joblib'
    
    joblib.dump(best_svm, model_filename)
    print(f"Successfully saved the trained SVM model to: {model_filename}")
    
    joblib.dump(scaler, scaler_filename)
    print(f"Successfully saved the feature scaler to: {scaler_filename}")

# --- Script Execution ---
if __name__ == "__main__":
    main()
