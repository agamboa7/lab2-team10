#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Final Model Comparison and Error Analysis Script

This script performs a comparative analysis between a pre-trained SVM model and
a Von Heijne model. It uses the training and prediction functions from the
'vonHeijne.py' script.

The analysis includes:
1.  Training a final Von Heijne model on the full training dataset.
2.  Loading a pre-trained SVM model and its feature scaler.
3.  Generating predictions from both models on the hold-out test set.
4.  Extending the test dataframe with these predictions.
5.  Calculating and comparing the False Positive Rate (FPR) on transmembrane proteins.
6.  Analyzing the characteristics of False Negative and False Positive predictions.
7.  Saving the combined results to a TSV file.

Prerequisites:
- A file named 'vonHeijne.py' must be in the same directory.
- 'merged_dataset_with_seqs.tsv', 'svm_signal_peptide_model.joblib', and
  'feature_scaler.joblib' must also be in the same directory.
"""
import pandas as pd
import numpy as np
import joblib
from sklearn.metrics import confusion_matrix
from collections import Counter

# --- Import functions from your existing script ---
# Ensure 'model_creation.py' is in the same directory
try:
    from vonHeijne import training as vh_training, prediction as vh_prediction
except ImportError:
    print("Error: 'vonHeijne.py' not found.")
    print("Please ensure your script containing the Von Heijne model is in the same directory.")
    exit()

# --- Constants and Functions for SVM Feature Engineering ---
# These are needed to process the test data for the SVM model
AMINO_ACIDS = sorted(list("ACDEFGHIKLMNPQRSTVWY"))
N_TERMINUS_LEN = 22
HYDROPHOBICITY_SCALE_KD = {
    'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5, 'Q': -3.5, 'E': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5, 'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6, 'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2
}
CHARGE_SCALE = {'D': -1, 'E': -1, 'K': 1, 'R': 1}
# Additional scales from the first script to ensure full feature set is created
TM_HELIX_PROPENSITY_SCALE = {
    'A': 0.25, 'R': -1.80, 'N': -0.64, 'D': -0.72, 'C': 0.04, 'Q': -0.69, 'E': -0.62, 'G': 0.16, 'H': -0.40, 'I': 0.73, 'L': 0.53, 'K': -1.10, 'M': 0.26, 'F': 0.61, 'P': -0.07, 'S': -0.26, 'T': -0.18, 'W': 0.37, 'Y': 0.02, 'V': 0.54
}
ALPHA_HELIX_PROPENSITY = {
    'A': 1.42, 'R': 0.98, 'N': 0.67, 'D': 1.01, 'C': 0.70, 'Q': 1.11, 'E': 1.51, 'G': 0.57, 'H': 1.00, 'I': 1.08, 'L': 1.21, 'K': 1.16, 'M': 1.45, 'F': 1.13, 'P': 0.57, 'S': 0.77, 'T': 0.83, 'W': 1.08, 'Y': 0.69, 'V': 1.06
}

def calculate_average_property(sequence, scale):
    if not sequence: return 0.0
    return sum(scale.get(aa, 0.0) for aa in sequence) / len(sequence)

def create_augmented_features(sequence):
    """Re-creates the exact feature set used for the SVM."""
    if not isinstance(sequence, str): sequence = ""
    n_term_seq = sequence[:N_TERMINUS_LEN].upper()
    n_region, h_region = n_term_seq[0:5], n_term_seq[5:16]
    if not n_term_seq: return [0.0] * (len(AMINO_ACIDS) + 6)
    counts, total_len = Counter(n_term_seq), len(n_term_seq)
    aa_composition = [(counts.get(aa, 0) / total_len) for aa in AMINO_ACIDS]
    n_region_charge = sum(CHARGE_SCALE.get(aa, 0) for aa in n_region)
    h_region_hydrophobicity_kd = calculate_average_property(h_region, HYDROPHOBICITY_SCALE_KD)
    h_region_propensity_tm = calculate_average_property(h_region, TM_HELIX_PROPENSITY_SCALE)
    avg_hydrophobicity_total = calculate_average_property(n_term_seq, HYDROPHOBICITY_SCALE_KD)
    avg_propensity_tm_total = calculate_average_property(n_term_seq, TM_HELIX_PROPENSITY_SCALE)
    avg_propensity_alpha_total = calculate_average_property(n_term_seq, ALPHA_HELIX_PROPENSITY)
    return aa_composition + [n_region_charge, h_region_hydrophobicity_kd, h_region_propensity_tm, avg_hydrophobicity_total, avg_propensity_tm_total, avg_propensity_alpha_total]

# --- Analysis Functions ---

def calculate_tm_fpr(df, model_name):
    """Calculates standard FPR and the specific FPR for the TM subset."""
    y_true = df['label']
    y_pred = df[f'predicted_{model_name.lower()}']

    cm = confusion_matrix(y_true, y_pred)
    # Handle cases where a class is not present in the predictions
    if cm.shape == (1, 1):
        if y_true.iloc[0] == 0: tn, fp, fn, tp = cm[0,0], 0, 0, 0
        else: tn, fp, fn, tp = 0, 0, 0, cm[0,0]
    else:
        tn, fp, fn, tp = cm.ravel()

    standard_fpr = fp / (fp + tn) if (fp + tn) > 0 else 0
    print(f"\n--- {model_name} Model ---")
    print(f"Standard FPR: {standard_fpr:.4f}")

    is_tm = df['Transmembrane_Helix_N_Terminus'] == True
    is_true_negative = df['label'] == 0
    negative_tm_df = df[is_tm & is_true_negative]

    if not negative_tm_df.empty:
        fp_tm = (df.loc[negative_tm_df.index, f'predicted_{model_name.lower()}'] == 1).sum()
        negative_tm_count = len(negative_tm_df)
        fpr_tm = fp_tm / negative_tm_count if negative_tm_count > 0 else 0
        print(f"Total Negative TM Proteins: {negative_tm_count}")
        print(f"False Positives from TM subset: {fp_tm}")
        print(f"FPR for Transmembrane proteins: {fpr_tm:.4f}")
    else:
        print("No transmembrane proteins found among the true negatives.")

def analyze_misclassifications(df, model_name):
    """Analyzes average features of misclassified sequences."""
    pred_col = f'predicted_{model_name.lower()}'
    is_fp = (df['label'] == 0) & (df[pred_col] == 1)
    is_fn = (df['label'] == 1) & (df[pred_col] == 0)
    is_tp = (df['label'] == 1) & (df[pred_col] == 1)
    is_tn = (df['label'] == 0) & (df[pred_col] == 0)

    def get_avg_features(df_subset):
        if df_subset.empty:
            return {'n_charge': 0, 'h_hydrophobicity': 0, 'count': 0}
        sequences = df_subset['Sequence']
        features = sequences.apply(lambda seq: (
            sum(CHARGE_SCALE.get(aa, 0) for aa in str(seq)[:5]),
            calculate_average_property(str(seq)[5:16], HYDROPHOBICITY_SCALE_KD)
        )).tolist()
        avg_charge = np.mean([f[0] for f in features])
        avg_hydro = np.mean([f[1] for f in features])
        return {'n_charge': avg_charge, 'h_hydrophobicity': avg_hydro, 'count': len(df_subset)}

    analysis_summary = pd.DataFrame({
        'True Positives': get_avg_features(df[is_tp]),
        'False Negatives': get_avg_features(df[is_fn]),
        'True Negatives': get_avg_features(df[is_tn]),
        'False Positives': get_avg_features(df[is_fp])
    }).T.round(3)

    print(f"\n--- {model_name} Model Misclassification Analysis ---")
    print("Average of key biological features by prediction outcome:")
    print(analysis_summary)
    print("\nInterpretation:")
    print("-> Compare FN hydrophobicity to TP: Lower values suggest the model missed less canonical signal peptides.")
    print("-> Compare FP hydrophobicity to TN: Higher values suggest the model was confused by hydrophobic N-termini (like TM anchors).")


def main():
    """Main function to run the entire comparison pipeline."""
    print("--- Starting Final Model Comparison and Error Analysis ---")

    # --- 1. Load Data ---
    try:
        df = pd.read_csv("merged_dataset_with_seqs.tsv", sep="\t")
        train_df = df[df["Type"] == "train"].copy()
        test_df = df[df["Type"] == "test"].copy()
        print(f"Data loaded. Training set: {len(train_df)}, Test set: {len(test_df)}")
    except FileNotFoundError:
        print("Error: 'merged_dataset_with_seqs.tsv' not found.")
        return

    # --- 2. Train and Predict with Final Von Heijne Model ---
    print("\n[1] Training final Von Heijne model on the full training set...")
    # Use the full training set for both training and validation to get a single, robust PSWM and threshold
    final_pswm, final_threshold, _, _ = vh_training(train_df, train_df)
    print(f"Final Von Heijne PSWM trained. Optimal threshold: {final_threshold:.4f}")
    print("Making predictions with Von Heijne model on the test set...")
    y_pred_vh = vh_prediction(test_df, final_pswm, final_threshold)

    # --- 3. Load and Predict with Final SVM Model ---
    print("\n[2] Loading pre-trained SVM model and making predictions...")
    try:
        svm_model = joblib.load('svm_signal_peptide_model.joblib')
        scaler = joblib.load('feature_scaler.joblib')
    except FileNotFoundError:
        print("Error: Could not find 'svm_signal_peptide_model.joblib' or 'feature_scaler.joblib'.")
        print("Please run the SVM training script first to generate these files.")
        return

    # Create features for the test set
    X_test_list = test_df['Sequence'].apply(create_augmented_features).tolist()
    X_test = np.array(X_test_list)
    X_test_scaled = scaler.transform(X_test)

    # The SVM was trained on the top 15 features. We must select the *same* 15 features for prediction.
    # We can find these indices by re-running feature selection on the training data.
    from sklearn.ensemble import RandomForestClassifier
    print("Identifying the top 15 features used by the SVM...")
    X_train_list = train_df['Sequence'].apply(create_augmented_features).tolist()
    X_train_scaled = scaler.transform(np.array(X_train_list))
    rf_selector = RandomForestClassifier(n_estimators=150, random_state=42, n_jobs=-1).fit(X_train_scaled, train_df['label'])
    indices = np.argsort(rf_selector.feature_importances_)[::-1][:15]
    
    X_test_selected = X_test_scaled[:, indices]
    print(f"Filtered test features to shape: {X_test_selected.shape}")
    y_pred_svm = svm_model.predict(X_test_selected)

    # --- 4. Extend DataFrame with Predictions ---
    print("\n[3] Extending test DataFrame with model predictions...")
    test_df_extended = test_df.copy()
    test_df_extended['predicted_svm'] = y_pred_svm
    test_df_extended['predicted_von_heijne'] = y_pred_vh
    output_filename = 'test_set_with_predictions.tsv'
    test_df_extended.to_csv(output_filename, sep='\t', index=False)
    print(f"Saved test set with predictions to '{output_filename}'")

    # --- 5. FPR Analysis for Transmembrane Proteins ---
    print("\n[4] Calculating FPR for Transmembrane (TM) Proteins...")
    calculate_tm_fpr(test_df_extended, 'SVM')
    calculate_tm_fpr(test_df_extended, 'Von_Heijne')
    print("\nComparison: A higher FPR for TM proteins suggests a model struggles to distinguish")
    print("their hydrophobic N-terminal anchors from true signal peptides.")

    # --- 6. Analysis of Misclassified Sequences ---
    print("\n[5] Analyzing False Negative and False Positive Predictions...")
    analyze_misclassifications(test_df_extended, 'SVM')
    analyze_misclassifications(test_df_extended, 'Von_Heijne')

    print("\n--- Analysis Complete ---")

if __name__ == '__main__':
    main()