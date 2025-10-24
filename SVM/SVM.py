#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
LB2 Project: Support Vector Machine for Signal Peptide Detection
Practical Session IV

This script implements a complete machine learning pipeline to classify proteins
based on the presence of a signal peptide (SP).

The pipeline includes:
1. Loading data with predefined cross-validation folds.
2. Engineering a rich set of features based on the physicochemical properties
   of the N-terminal sequence.
3. Applying feature selection using a RandomForestClassifier to rank features.
4. Filtering the dataset to keep only the top 15 most important features.
5. Performing hyperparameter optimization for an SVM classifier using GridSearchCV,
   with a feature to save and load the best parameters to avoid re-running.
6. Evaluating the final, optimized model on a hold-out test set.
7. Saving the trained model and the feature scaler to disk for future use.

Usage:
    python your_script_name.py --input merged_dataset_with_seqs.tsv

"""

import argparse
import json
import os
from collections import Counter

import joblib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import (auc, classification_report, confusion_matrix,
                             roc_curve)
from sklearn.model_selection import GridSearchCV
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC

# --- Constants and Configuration ---
# Set plot style for consistency
sns.set_theme(context="notebook", style="darkgrid", palette="muted")

# Feature Engineering Constants
AMINO_ACIDS = sorted(list("ACDEFGHIKLMNPQRSTVWY"))
N_TERMINUS_LEN = 22
HYDROPHOBICITY_SCALE_KD = {
    'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5, 'Q': -3.5, 'E': -3.5,
    'G': -0.4, 'H': -3.2, 'I': 4.5, 'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8,
    'P': -1.6, 'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2
}
TM_HELIX_PROPENSITY_SCALE = {
    'A': 0.25, 'R': -1.80, 'N': -0.64, 'D': -0.72, 'C': 0.04, 'Q': -0.69,
    'E': -0.62, 'G': 0.16, 'H': -0.40, 'I': 0.73, 'L': 0.53, 'K': -1.10,
    'M': 0.26, 'F': 0.61, 'P': -0.07, 'S': -0.26, 'T': -0.18, 'W': 0.37,
    'Y': 0.02, 'V': 0.54
}
ALPHA_HELIX_PROPENSITY = {
    'A': 1.42, 'R': 0.98, 'N': 0.67, 'D': 1.01, 'C': 0.70, 'Q': 1.11,
    'E': 1.51, 'G': 0.57, 'H': 1.00, 'I': 1.08, 'L': 1.21, 'K': 1.16,
    'M': 1.45, 'F': 1.13, 'P': 0.57, 'S': 0.77, 'T': 0.83, 'W': 1.08,
    'Y': 0.69, 'V': 1.06
}
CHARGE_SCALE = {'D': -1, 'E': -1, 'K': 1, 'R': 1}

# File paths
BEST_PARAMS_FILENAME = 'best_svm_params.json'
MODEL_FILENAME = 'svm_signal_peptide_model.joblib'
SCALER_FILENAME = 'feature_scaler.joblib'


def load_data(filepath):
    """Loads the dataset and splits it into training and testing sets."""
    print("--- 1. Loading and Preparing Data ---")
    try:
        df = pd.read_csv(filepath, sep="\t")
        print(f"Dataset '{filepath}' loaded successfully.")
    except FileNotFoundError:
        print(f"Error: Dataset file not found at '{filepath}'.")
        exit()

    train_df = df[df["Type"] == "train"].copy()
    test_df = df[df["Type"] == "test"].copy()
    print(f"Training set size: {len(train_df)}")
    print(f"Test set size: {len(test_df)}")
    return train_df, test_df


def create_augmented_features(sequence):
    """Engineers a feature vector from a single protein sequence."""
    if not isinstance(sequence, str):
        sequence = ""

    n_term_seq = sequence[:N_TERMINUS_LEN].upper()
    n_region, h_region = n_term_seq[0:5], n_term_seq[5:16]

    if not n_term_seq:
        # Return a zero vector of the correct length if sequence is empty
        return [0.0] * (len(AMINO_ACIDS) + 6)

    # Helper for calculating average properties
    def calculate_average_property(seq, scale):
        if not seq:
            return 0.0
        return sum(scale.get(aa, 0.0) for aa in seq) / len(seq)

    # Feature calculations
    counts, total_len = Counter(n_term_seq), len(n_term_seq)
    aa_composition = [(counts.get(aa, 0) / total_len) for aa in AMINO_ACIDS]
    n_region_charge = sum(CHARGE_SCALE.get(aa, 0) for aa in n_region)
    h_region_hydrophobicity_kd = calculate_average_property(h_region, HYDROPHOBICITY_SCALE_KD)
    h_region_propensity_tm = calculate_average_property(h_region, TM_HELIX_PROPENSITY_SCALE)
    avg_hydrophobicity_total = calculate_average_property(n_term_seq, HYDROPHOBICITY_SCALE_KD)
    avg_propensity_tm_total = calculate_average_property(n_term_seq, TM_HELIX_PROPENSITY_SCALE)
    avg_propensity_alpha_total = calculate_average_property(n_term_seq, ALPHA_HELIX_PROPENSITY)

    return (aa_composition + [
        n_region_charge, h_region_hydrophobicity_kd, h_region_propensity_tm,
        avg_hydrophobicity_total, avg_propensity_tm_total, avg_propensity_alpha_total
    ])


def process_features(df, scaler=None):
    """Applies feature engineering and scaling to a dataframe."""
    print(f"\n--- 2. Engineering Features for {'Training' if scaler is None else 'Test'} Data ---")
    
    X_list = df['Sequence'].apply(create_augmented_features).tolist()
    X = np.array(X_list)
    y = df['label'].values

    feature_names = ([f"AA_{aa}_freq" for aa in AMINO_ACIDS] + [
        "n_region_charge", "h_region_hydrophobicity", "h_region_TM_propensity",
        "total_avg_hydrophobicity", "total_avg_TM_propensity", "total_avg_alpha_propensity"
    ])

    if scaler is None:
        scaler = StandardScaler()
        X_scaled = scaler.fit_transform(X)
        print(f"Training feature engineering complete. Shape: {X_scaled.shape}")
        return X_scaled, y, feature_names, scaler
    else:
        X_scaled = scaler.transform(X)
        print(f"Test feature engineering complete. Shape: {X_scaled.shape}")
        return X_scaled, y


def select_features(X_train_scaled, y_train, feature_names):
    """Selects the most important features using a RandomForest."""
    print("\n--- 3. Running Feature Selection with RandomForest ---")
    rf_selector = RandomForestClassifier(n_estimators=150, random_state=42, n_jobs=-1)
    rf_selector.fit(X_train_scaled, y_train)

    importances = rf_selector.feature_importances_
    feature_importance_df = pd.DataFrame({
        'Feature': feature_names,
        'Importance': importances
    }).sort_values(by='Importance', ascending=False)

    # Plot and save feature importances
    plt.figure(figsize=(12, 10))
    sns.barplot(x='Importance', y='Feature', data=feature_importance_df,
                hue='Feature', palette='viridis_r', legend=False)
    plt.title('Augmented Feature Importances from RandomForest', fontsize=16)
    plt.tight_layout()
    plt.savefig("svm_augmented_feature_importances.png", dpi=300)
    plt.close() # Close plot to prevent display in script mode

    print("Top 15 most important features:\n", feature_importance_df.head(15).to_string())

    top_15_feature_names = feature_importance_df.head(15)['Feature'].tolist()
    top_15_indices = [feature_names.index(feature) for feature in top_15_feature_names]

    return top_15_indices


def train_svm(X_train_selected, y_train, train_df):
    """Optimizes SVM hyperparameters using GridSearchCV with caching."""
    print("\n--- 4. Optimizing SVM Hyperparameters with GridSearchCV ---")

    if os.path.exists(BEST_PARAMS_FILENAME):
        print(f"Found existing parameters in '{BEST_PARAMS_FILENAME}'. Loading.")
        with open(BEST_PARAMS_FILENAME, 'r') as f:
            best_params = json.load(f)
        print(f"Loaded parameters: {best_params}")
        
        best_svm = SVC(**best_params, probability=True, random_state=42)
        best_svm.fit(X_train_selected, y_train)
        print("Model training with pre-loaded parameters complete.")

    else:
        print(f"'{BEST_PARAMS_FILENAME}' not found. Running GridSearchCV.")
        param_grid_nested = [
            {'kernel': ['linear'], 'C': [0.1, 1, 10, 100]},
            {'kernel': ['rbf'], 'C': [0.1, 1, 10, 100], 'gamma': ['scale', 'auto', 0.01, 0.1]},
            {'kernel': ['poly'], 'C': [0.1, 1, 10, 100], 'degree': [2, 3], 'gamma': ['scale', 'auto']}
        ]

        # Create custom cross-validation iterator from folds
        cv_iterator = []
        index_mapper = {original_idx: i for i, original_idx in enumerate(train_df.index)}
        for i in range(5):
            train_fold_indices = train_df[train_df['fold'] != i].index
            val_fold_indices = train_df[train_df['fold'] == i].index
            train_indices_mapped = [index_mapper[idx] for idx in train_fold_indices]
            val_indices_mapped = [index_mapper[idx] for idx in val_fold_indices]
            cv_iterator.append((train_indices_mapped, val_indices_mapped))

        grid_search = GridSearchCV(
            estimator=SVC(probability=True, random_state=42),
            param_grid=param_grid_nested,
            cv=cv_iterator,
            scoring='f1',
            n_jobs=-1,
            verbose=2
        )
        grid_search.fit(X_train_selected, y_train)
        
        best_svm = grid_search.best_estimator_
        print("\nGridSearchCV complete.")
        print(f"Best parameters found: {grid_search.best_params_}")
        print(f"Best cross-validation F1-score: {grid_search.best_score_:.4f}")
        
        print(f"\n--- Saving Best Parameters to {BEST_PARAMS_FILENAME} ---")
        with open(BEST_PARAMS_FILENAME, 'w') as f:
            json.dump(grid_search.best_params_, f, indent=4)
        print("Successfully saved the best parameters.")
        
    return best_svm


def evaluate_model(model, X_test_selected, y_test):
    """Evaluates the final model on the test set and saves plots."""
    print("\n--- 5. Evaluating the Best Model on the Hold-out Test Set ---")

    if X_test_selected.shape[0] == 0:
        print("\nTest set is empty. Skipping final evaluation.")
        return

    y_pred = model.predict(X_test_selected)
    y_pred_proba = model.predict_proba(X_test_selected)[:, 1]

    print("\nClassification Report on Test Set:")
    print(classification_report(y_test, y_pred, target_names=['Negative (0)', 'Positive (1)']))

    # Confusion Matrix
    cm = confusion_matrix(y_test, y_pred)
    plt.figure(figsize=(8, 6))
    sns.heatmap(cm, annot=True, fmt='d', cmap='Blues', cbar=False,
                xticklabels=['Negative', 'Positive'], yticklabels=['Negative', 'Positive'])
    plt.title('Confusion Matrix on Test Set', fontsize=16)
    plt.xlabel('Predicted Label', fontsize=12)
    plt.ylabel('True Label', fontsize=12)
    plt.savefig("svm_confusion_matrix.png", dpi=300)
    plt.close()

    # ROC Curve
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
    plt.savefig("svm_roc_curve.png", dpi=300)
    plt.close()

    print("\nEvaluation complete. Plots saved to 'svm_confusion_matrix.png' and 'svm_roc_curve.png'.")


def main(args):
    """Main function to run the complete ML pipeline."""
    # 1. Load data
    train_df, test_df = load_data(args.input)

    # 2. Engineer features for training data
    X_train_scaled, y_train, feature_names, scaler = process_features(train_df)

    # 3. Select top features
    top_15_indices = select_features(X_train_scaled, y_train, feature_names)
    
    # 3.5 Filter dataset to top features
    print("\n--- 3.5. Filtering Dataset to Top 15 Features ---")
    X_train_selected = X_train_scaled[:, top_15_indices]
    print(f"Filtered training features shape: {X_train_selected.shape}")

    # 4. Train the SVM model
    best_svm = train_svm(X_train_selected, y_train, train_df)

    # 5. Process test data and evaluate
    if not test_df.empty:
        X_test_scaled, y_test = process_features(test_df, scaler=scaler)
        X_test_selected = X_test_scaled[:, top_15_indices]
        print(f"Filtered test features shape: {X_test_selected.shape}")
        evaluate_model(best_svm, X_test_selected, y_test)
    else:
        print("\nNo test data to evaluate.")

    # 6. Save final model and scaler
    print("\n--- 6. Saving the final trained model and scaler to disk ---")
    joblib.dump(best_svm, MODEL_FILENAME)
    print(f"Successfully saved the trained SVM model to: {MODEL_FILENAME}")
    joblib.dump(scaler, SCALER_FILENAME)
    print(f"Successfully saved the feature scaler to: {SCALER_FILENAME}")

    print("\n--- Pipeline Finished ---")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="SVM Signal Peptide Detection Pipeline.")
    parser.add_argument(
        "--input",
        type=str,
        required=True,
        help="Path to the input 'merged_dataset_with_seqs.tsv' file."
    )
    args = parser.parse_args()
    main(args)