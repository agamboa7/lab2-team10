#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Final Model Comparison: SVM vs. Deep Learning vs. Von Heijne

This script performs a head-to-head comparative analysis between three pre-trained models:
1.  Support Vector Machine (SVM)
2.  PyTorch Deep Learning (CNN + BiLSTM)
3.  Von Heijne Position-Specific Weight Matrix (PSWM)

The analysis includes:
1.  Loading the hold-out test set.
2.  Loading the pre-trained SVM model and its feature scaler.
3.  Loading the pre-trained Deep Learning model architecture and weights.
4.  Training a final Von Heijne model on the full training dataset.
5.  Generating predictions from all three models on the test set.
6.  Comparing overall performance metrics (Classification Report) in a markdown file.
7.  Comparing the False Positive Rate (FPR) on transmembrane proteins for each model.
8.  Analyzing the characteristics of misclassified predictions for each model.
9.  Saving the combined results to a TSV file.

Prerequisites:
- 'merged_dataset_with_seqs.tsv'
- 'svm_signal_peptide_model.joblib' (from SVM.py)
- 'feature_scaler.joblib' (from SVM.py)
- 'best_signal_peptide_model.pth' (from your DL script)
- 'vonHeijne.py' (containing the Von Heijne training and prediction functions)
"""
import pandas as pd
import numpy as np
import joblib
import torch
import torch.nn as nn
from torch.utils.data import Dataset, DataLoader
from sklearn.metrics import confusion_matrix, classification_report
from collections import Counter
import warnings

# --- Import Von Heijne functions ---
try:
    from vonHeijne import training as vh_training, prediction as vh_prediction
except ImportError:
    print("Error: 'vonHeijne.py' not found.")
    print("Please ensure the script containing the Von Heijne model is in the same directory.")
    exit()

warnings.filterwarnings('ignore', category=UserWarning, module='sklearn')

# =============================================================================
# --- SECTION 1: DEEP LEARNING MODEL DEFINITIONS ---
# =============================================================================

class DLConfig:
    """Minimal config needed for the DL model architecture."""
    INPUT_SIZE = 20
    CONV_OUT_CHANNELS = 64
    CONV_KERNEL_SIZE = 17
    LSTM_HIDDEN_SIZE = 256
    LSTM_NUM_LAYERS = 2
    DROPOUT_PROB = 0.5
    MAX_SEQ_LEN = 70
    BATCH_SIZE = 32
    DEVICE = "cuda" if torch.cuda.is_available() else "cpu"

class SignalPeptideClassifier(nn.Module):
    """CNN + BiLSTM model for sequence classification."""
    def __init__(self, config: DLConfig):
        super().__init__()
        self.conv_block = nn.Sequential(
            nn.Conv1d(in_channels=config.INPUT_SIZE, out_channels=config.CONV_OUT_CHANNELS, kernel_size=config.CONV_KERNEL_SIZE, padding='same'),
            nn.BatchNorm1d(config.CONV_OUT_CHANNELS), nn.ReLU()
        )
        self.lstm = nn.LSTM(
            input_size=config.CONV_OUT_CHANNELS, hidden_size=config.LSTM_HIDDEN_SIZE, num_layers=config.LSTM_NUM_LAYERS,
            batch_first=True, bidirectional=True, dropout=config.DROPOUT_PROB if config.LSTM_NUM_LAYERS > 1 else 0
        )
        self.fc_block = nn.Sequential(
            nn.Linear(config.LSTM_HIDDEN_SIZE * 2, 512), nn.ReLU(), nn.Dropout(config.DROPOUT_PROB),
            nn.Linear(512, 256), nn.ReLU(), nn.Dropout(config.DROPOUT_PROB),
            nn.Linear(256, 1)
        )
    def forward(self, x):
        x = x.permute(0, 2, 1); x = self.conv_block(x); x = x.permute(0, 2, 1)
        out, _ = self.lstm(x); out = out[:, -1, :]
        return self.fc_block(out).squeeze(-1)

def one_hot_encode_dl(seq: str, max_len: int, aa_map: dict) -> np.ndarray:
    one_hot = np.zeros((max_len, len(aa_map)), dtype=np.float32)
    for i, aa in enumerate(seq[:max_len]):
        if aa in aa_map: one_hot[i, aa_map[aa]] = 1.0
    return one_hot

class DLDataset(Dataset):
    def __init__(self, sequences): self.sequences = sequences
    def __len__(self): return len(self.sequences)
    def __getitem__(self, idx): return torch.tensor(self.sequences[idx], dtype=torch.float32)

def get_dl_predictions(sequences, model_path, config):
    """Loads a DL model and generates predictions for a list of sequences."""
    print("Making predictions with Deep Learning model...")
    aa_alph = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    aa_map = {aa: i for i, aa in enumerate(aa_alph)}

    model = SignalPeptideClassifier(config).to(config.DEVICE)
    model.load_state_dict(torch.load(model_path, map_location=torch.device(config.DEVICE)))
    model.eval()

    encoded_seqs = [one_hot_encode_dl(s, config.MAX_SEQ_LEN, aa_map) for s in sequences]
    dataset = DLDataset(encoded_seqs)
    loader = DataLoader(dataset, batch_size=config.BATCH_SIZE, shuffle=False)

    all_preds = []
    with torch.no_grad():
        for batch_seqs in loader:
            batch_seqs = batch_seqs.to(config.DEVICE)
            logits = model(batch_seqs)
            preds = (torch.sigmoid(logits) > 0.5).int()
            all_preds.extend(preds.cpu().numpy())
    return np.array(all_preds)

# =============================================================================
# --- SECTION 2: SVM FEATURE ENGINEERING DEFINITIONS ---
# =============================================================================

AMINO_ACIDS = sorted(list("ACDEFGHIKLMNPQRSTVWY"))
N_TERMINUS_LEN = 22
HYDROPHOBICITY_SCALE_KD = {'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5, 'Q': -3.5, 'E': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5, 'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6, 'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2}
TM_HELIX_PROPENSITY_SCALE = {'A': 0.25, 'R': -1.80, 'N': -0.64, 'D': -0.72, 'C': 0.04, 'Q': -0.69, 'E': -0.62, 'G': 0.16, 'H': -0.40, 'I': 0.73, 'L': 0.53, 'K': -1.10, 'M': 0.26, 'F': 0.61, 'P': -0.07, 'S': -0.26, 'T': -0.18, 'W': 0.37, 'Y': 0.02, 'V': 0.54}
ALPHA_HELIX_PROPENSITY = {'A': 1.42, 'R': 0.98, 'N': 0.67, 'D': 1.01, 'C': 0.70, 'Q': 1.11, 'E': 1.51, 'G': 0.57, 'H': 1.00, 'I': 1.08, 'L': 1.21, 'K': 1.16, 'M': 1.45, 'F': 1.13, 'P': 0.57, 'S': 0.77, 'T': 0.83, 'W': 1.08, 'Y': 0.69, 'V': 1.06}
CHARGE_SCALE = {'D': -1, 'E': -1, 'K': 1, 'R': 1}

def calculate_average_property(sequence, scale):
    if not sequence: return 0.0
    return sum(scale.get(aa, 0.0) for aa in sequence) / len(sequence)

def create_augmented_features(sequence):
    if not isinstance(sequence, str): sequence = ""
    n_term_seq = sequence[:N_TERMINUS_LEN].upper()
    n_region, h_region = n_term_seq[0:5], n_term_seq[5:16]
    if not n_term_seq: return [0.0] * (len(AMINO_ACIDS) + 6)
    counts, total_len = Counter(n_term_seq), len(n_term_seq)
    aa_composition = [(counts.get(aa, 0) / total_len) for aa in AMINO_ACIDS]
    n_charge = sum(CHARGE_SCALE.get(aa, 0) for aa in n_region)
    h_hydro = calculate_average_property(h_region, HYDROPHOBICITY_SCALE_KD)
    h_tm = calculate_average_property(h_region, TM_HELIX_PROPENSITY_SCALE)
    avg_hydro = calculate_average_property(n_term_seq, HYDROPHOBICITY_SCALE_KD)
    avg_tm = calculate_average_property(n_term_seq, TM_HELIX_PROPENSITY_SCALE)
    avg_alpha = calculate_average_property(n_term_seq, ALPHA_HELIX_PROPENSITY)
    return aa_composition + [n_charge, h_hydro, h_tm, avg_hydro, avg_tm, avg_alpha]

# =============================================================================
# --- SECTION 3: COMPARATIVE ANALYSIS FUNCTIONS ---
# =============================================================================

def classification_report_to_markdown(report_dict):
    """Converts a scikit-learn classification report dictionary to a Markdown table."""
    header = "| Metric          | Precision | Recall | F1-Score | Support |"
    separator = "|-----------------|-----------|--------|----------|---------|"
    
    lines = [header, separator]
    
    for label, metrics in report_dict.items():
        if isinstance(metrics, dict) and label in ['Negative', 'Positive']:
            row = f"| **{label}**       | {metrics['precision']:.2f}      | {metrics['recall']:.2f}   | {metrics['f1-score']:.2f}     | {metrics['support']}      |"
            lines.append(row)
    
    lines.append(separator)
    acc = report_dict.get('accuracy', 0)
    lines.append(f"| **Accuracy**      |           |        | {acc:.2f}     | {report_dict.get('macro avg', {}).get('support')}      |")
    
    macro_avg = report_dict.get('macro avg', {})
    lines.append(f"| **Macro Avg**     | {macro_avg.get('precision', 0):.2f}      | {macro_avg.get('recall', 0):.2f}   | {macro_avg.get('f1-score', 0):.2f}     | {macro_avg.get('support')}      |")
    
    weighted_avg = report_dict.get('weighted avg', {})
    lines.append(f"| **Weighted Avg**  | {weighted_avg.get('precision', 0):.2f}      | {weighted_avg.get('recall', 0):.2f}   | {weighted_avg.get('f1-score', 0):.2f}     | {weighted_avg.get('support')}      |")
    
    return "\n".join(lines)

def calculate_tm_fpr(df, model_name, file=None):
    """Calculates the specific FPR for the Transmembrane protein subset."""
    pred_col = f'predicted_{model_name.lower()}'
    
    def write_line(text):
        print(text)
        if file: file.write(text + '\n')

    write_line(f"\n### {model_name} Model: Transmembrane Protein Analysis")
    
    # Find the correct transmembrane column name
    possible_cols = ['Transmembrane_Helix_N_Terminus', 'Transmembrane', 'Transmem']
    tm_col = next((c for c in possible_cols if c in df.columns), None)
    if tm_col is None:
        matches = [c for c in df.columns if 'transmem' in c.lower()]
        tm_col = matches[0] if matches else None
    
    if tm_col is None:
        raise KeyError(f"No transmembrane column found. Available columns: {list(df.columns)}")

    # Ensure boolean conversion is robust
    col_vals = df[tm_col]
    if pd.api.types.is_bool_dtype(col_vals):
        is_tm = col_vals
    else:
        is_tm = col_vals.astype(str).str.lower().isin(['true', '1', 't', 'yes'])

    is_true_negative = df['label'] == 0
    negative_tm_df = df[is_tm & is_true_negative]

    if not negative_tm_df.empty:
        fp_tm = (df.loc[negative_tm_df.index, pred_col] == 1).sum()
        negative_tm_count = len(negative_tm_df)
        fpr_tm = fp_tm / negative_tm_count if negative_tm_count > 0 else 0
        write_line(f"- **Total Negative TM Proteins:** {negative_tm_count}")
        write_line(f"- **False Positives from TM subset:** {fp_tm}")
        write_line(f"- **FPR for Transmembrane proteins:** `{fpr_tm:.4f}`")
    else:
        write_line("- No transmembrane proteins found among the true negatives.")

def analyze_misclassifications(df, model_name, file=None):
    """Analyzes average features of misclassified sequences."""
    pred_col = f'predicted_{model_name.lower()}'
    is_fp = (df['label'] == 0) & (df[pred_col] == 1)
    is_fn = (df['label'] == 1) & (df[pred_col] == 0)
    is_tp = (df['label'] == 1) & (df[pred_col] == 1)
    
    def get_avg_hydrophobicity(df_subset):
        if df_subset.empty: return 0
        return df_subset['Sequence'].apply(lambda s: calculate_average_property(str(s)[5:16], HYDROPHOBICITY_SCALE_KD)).mean()

    def write_line(text):
        print(text)
        if file: file.write(text + '\n')

    write_line(f"\n### {model_name} Model: Misclassification Analysis")
    tp_hydro = get_avg_hydrophobicity(df[is_tp])
    fn_hydro = get_avg_hydrophobicity(df[is_fn])
    fp_hydro = get_avg_hydrophobicity(df[is_fp])
    
    write_line(f"- **Average H-region Hydrophobicity of True Positives:** `{tp_hydro:.3f}`")
    write_line(f"- **Average H-region Hydrophobicity of FALSE NEGATIVES:** `{fn_hydro:.3f}`")
    write_line(f"- **Average H-region Hydrophobicity of FALSE POSITIVES:** `{fp_hydro:.3f}`")
    write_line("\n**Interpretation:**")
    if fn_hydro < tp_hydro:
        write_line("- FN hydrophobicity is lower than TP, suggesting the model misses less canonical (less hydrophobic) signal peptides.")
    write_line("- FP hydrophobicity > 0 suggests the model is confused by hydrophobic N-termini (like TM anchors).")

# =============================================================================
# --- SECTION 4: MAIN EXECUTION PIPELINE ---
# =============================================================================

def main():
    """Main function to run the entire comparison pipeline."""
    report_filename = "final_model_comparison_report.md"
    print(f"--- Starting Final Model Comparison: SVM vs. Deep Learning vs. Von Heijne ---")
    print(f"A detailed report will be saved to '{report_filename}'")

    with open(report_filename, "w") as report_file:
        def write_line(text, console_only=False):
            print(text.replace('`', '')) # Clean markdown for console
            if not console_only:
                report_file.write(text + '\n')

        write_line("# Final Model Comparison: SVM vs. Deep Learning vs. Von Heijne\n")

        # --- 1. Load Data ---
        write_line("## 1. Data Loading")
        try:
            df = pd.read_csv("merged_dataset_with_seqs.tsv", sep="\t")
            train_df = df[df["Type"] == "train"].copy()
            test_df = df[df["Type"] == "test"].copy()
            write_line(f"- Data loaded from `merged_dataset_with_seqs.tsv`.")
            write_line(f"- **Training set size:** {len(train_df)}")
            write_line(f"- **Test set size:** {len(test_df)}")
        except FileNotFoundError:
            write_line("Error: `merged_dataset_with_seqs.tsv` not found. Aborting.")
            return

        # --- 2. SVM Model Prediction ---
        write_line("\n## 2. Model Predictions")
        try:
            print("\n[SVM] Loading pre-trained model and making predictions...")
            svm_model = joblib.load('svm_signal_peptide_model.joblib')
            scaler = joblib.load('feature_scaler.joblib')
            
            X_test_svm_features = np.array(test_df['Sequence'].apply(create_augmented_features).tolist())
            X_test_scaled = scaler.transform(X_test_svm_features)

            # Re-run feature selection on training data to get the exact feature indices
            from sklearn.ensemble import RandomForestClassifier
            X_train_svm_features = np.array(train_df['Sequence'].apply(create_augmented_features).tolist())
            X_train_scaled = scaler.transform(X_train_svm_features)
            rf_selector = RandomForestClassifier(n_estimators=150, random_state=42, n_jobs=-1).fit(X_train_scaled, train_df['label'])
            top_15_indices = np.argsort(rf_selector.feature_importances_)[::-1][:15]
            
            X_test_selected = X_test_scaled[:, top_15_indices]
            y_pred_svm = svm_model.predict(X_test_selected)
            write_line("- SVM predictions generated.", console_only=True)
        except FileNotFoundError:
            write_line("Error: Could not find SVM model files. Skipping SVM analysis.")
            y_pred_svm = None

        # --- 3. Deep Learning Model Prediction ---
        try:
            print("\n[Deep Learning] Loading pre-trained model and making predictions...")
            dl_config = DLConfig()
            y_pred_dl = get_dl_predictions(
                sequences=test_df['Sequence'].tolist(),
                model_path='best_signal_peptide_model.pth',
                config=dl_config
            )
            write_line("- Deep Learning predictions generated.", console_only=True)
        except FileNotFoundError:
            write_line("Error: `best_signal_peptide_model.pth` not found. Skipping DL analysis.")
            y_pred_dl = None

        # --- 4. Von Heijne Model Prediction ---
        print("\n[Von Heijne] Training final model and making predictions...")
        final_pswm, final_threshold, _, _ = vh_training(train_df, train_df)
        y_pred_vh = vh_prediction(test_df, final_pswm, final_threshold)
        write_line(f"- Von Heijne model trained (threshold={final_threshold:.3f}) and predictions generated.", console_only=True)

        # --- 5. Combine and Save Results ---
        write_line("\n## 3. Saving Combined Results")
        test_df_extended = test_df.copy()
        if y_pred_svm is not None: test_df_extended['predicted_svm'] = y_pred_svm
        if y_pred_dl is not None: test_df_extended['predicted_dl'] = y_pred_dl
        test_df_extended['predicted_von_heijne'] = y_pred_vh
        
        output_filename = 'test_set_with_all_predictions.tsv'
        test_df_extended.to_csv(output_filename, sep='\t', index=False)
        write_line(f"- Extended test set with all model predictions saved to `{output_filename}`.")

        # --- 6. Comparative Analysis ---
        write_line("\n## 4. Comparative Performance Analysis")
        y_true = test_df_extended['label']
        
        models_to_analyze = []
        if 'predicted_svm' in test_df_extended: models_to_analyze.append(('SVM', y_pred_svm))
        if 'predicted_dl' in test_df_extended: models_to_analyze.append(('DL', y_pred_dl))
        models_to_analyze.append(('Von_Heijne', y_pred_vh))
        
        for name, preds in models_to_analyze:
            write_line(f"\n### {name} Model: Classification Report")
            report_dict = classification_report(y_true, preds, target_names=['Negative', 'Positive'], output_dict=True)
            report_md = classification_report_to_markdown(report_dict)
            write_line(report_md)

        write_line("\n## 5. Specialized Error Analysis")
        for name, _ in models_to_analyze:
            calculate_tm_fpr(test_df_extended, name, file=report_file)
        
        for name, _ in models_to_analyze:
            analyze_misclassifications(test_df_extended, name, file=report_file)

    print("\n--- Analysis Complete ---")

if __name__ == '__main__':
    main()