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
7.  Generating visualizations for in-depth error analysis for both models, including:
    - Taxonomy distribution of errors.
    - Sequence logos of True Positives vs. False Negatives (with improved aesthetics).
    - Sequence length distributions.
    - Amino acid composition comparisons.
    - Physicochemical property comparisons.
8.  Saving the combined results to a TSV file.

Prerequisites:
- A file named 'vonHeijne.py' must be in the same directory.
- 'merged_dataset_with_seqs.tsv', 'svm_signal_peptide_model.joblib', and
  'feature_scaler.joblib' must also be in the same directory.
- Required libraries: pandas, numpy, scikit-learn, joblib, matplotlib, seaborn, logomaker.
"""
import pandas as pd
import numpy as np
import joblib
from sklearn.metrics import confusion_matrix
from collections import Counter
import os

# --- Import visualization libraries ---
try:
    import matplotlib.pyplot as plt
    import seaborn as sns
    import logomaker as lm
except ImportError:
    print("Error: Visualization libraries not found.")
    print("Please install them using: pip install matplotlib seaborn logomaker")
    exit()


# --- Import functions from your existing script ---
# Ensure 'vonHeijne.py' is in the same directory
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

# --- Physicochemical properties for new plots ---
AA_PHYSICAL_PROPERTIES = {
    'Hydrophobic': ['A', 'V', 'I', 'L', 'M', 'F', 'Y', 'W'],
    'Polar': ['S', 'T', 'N', 'Q', 'C', 'G', 'P'],
    'Charged (+)': ['R', 'K', 'H'],
    'Charged (-)': ['D', 'E']
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


# --- VISUALIZATION FUNCTIONS ---

def plot_taxonomy_pie_charts(fp_df, fn_df, model_name, output_dir):
    """Generates and saves pie charts for the taxonomic distribution of FP and FN."""
    if fp_df.empty and fn_df.empty:
        return

    fig, axes = plt.subplots(1, 2, figsize=(16, 8))
    fig.suptitle(f'Taxonomic Distribution of Misclassifications for {model_name} Model', fontsize=16)

    # False Positives
    if not fp_df.empty:
        fp_counts = fp_df['Kingdom'].value_counts()
        axes[0].pie(fp_counts, labels=fp_counts.index, autopct='%1.1f%%', startangle=140, colors=sns.color_palette("pastel"))
        axes[0].set_title(f'False Positives (n={len(fp_df)})')
    else:
        axes[0].text(0.5, 0.5, 'No False Positives', ha='center', va='center')
        axes[0].set_title('False Positives')
        axes[0].axis('off')

    # False Negatives
    if not fn_df.empty:
        fn_counts = fn_df['Kingdom'].value_counts()
        axes[1].pie(fn_counts, labels=fn_counts.index, autopct='%1.1f%%', startangle=140, colors=sns.color_palette("pastel"))
        axes[1].set_title(f'False Negatives (n={len(fn_df)})')
    else:
        axes[1].text(0.5, 0.5, 'No False Negatives', ha='center', va='center')
        axes[1].set_title('False Negatives')
        axes[1].axis('off')

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    filename = os.path.join(output_dir, f'{model_name}_taxonomy_error_distribution.png')
    plt.savefig(filename)
    plt.close()
    print(f"Saved taxonomy pie charts to '{filename}'")

def plot_sequence_logos(tp_df, fn_df, model_name, output_dir):
    """Generates and saves sequence logos comparing TP and FN N-termini with y-axis headroom."""
    if tp_df.empty or fn_df.empty:
        print("Skipping logo generation as True Positives or False Negatives are empty.")
        return

    tp_seqs = [seq[:N_TERMINUS_LEN] for seq in tp_df['Sequence']]
    fn_seqs = [seq[:N_TERMINUS_LEN] for seq in fn_df['Sequence']]

    # Convert sequences to matrices
    tp_matrix = pd.DataFrame(lm.alignment_to_matrix(tp_seqs))
    fn_matrix = pd.DataFrame(lm.alignment_to_matrix(fn_seqs))

    # Create figure and axes
    fig, axes = plt.subplots(2, 1, figsize=(15, 6))

    # Draw the TP logo
    lm.Logo(tp_matrix, ax=axes[0], font_name='Arial Rounded MT Bold')
    axes[0].set_title(f'True Positives (n={len(tp_seqs)}) - {model_name}', fontsize=14)
    axes[0].set_ylabel('Bits')

    # Add y-axis headroom (e.g., +20%)
    ymin, ymax = axes[0].get_ylim()
    axes[0].set_ylim(ymin, ymax * 1.2)

    # Draw the FN logo
    lm.Logo(fn_matrix, ax=axes[1], font_name='Arial Rounded MT Bold')
    axes[1].set_title(f'False Negatives (n={len(fn_seqs)}) - {model_name}', fontsize=14)
    axes[1].set_xlabel('Position in N-Terminus')
    axes[1].set_ylabel('Bits')

    # Add y-axis headroom (e.g., +20%)
    ymin, ymax = axes[1].get_ylim()
    axes[1].set_ylim(ymin, ymax * 1.2)

    plt.tight_layout()
    filename = os.path.join(output_dir, f'{model_name}_TP_vs_FN_sequence_logos.png')
    plt.savefig(filename, dpi=300)
    plt.close()
    print(f"Saved sequence logos to '{filename}'")


def plot_length_comparison(tp_df, fn_df, model_name, output_dir):
    """Generates a box plot comparing sequence lengths of TPs and FNs."""
    if tp_df.empty or fn_df.empty:
        return

    tp_df['length'] = tp_df['Sequence'].str.len()
    fn_df['length'] = fn_df['Sequence'].str.len()

    plot_data = pd.concat([
        tp_df[['length']].assign(Category='True Positive'),
        fn_df[['length']].assign(Category='False Negative')
    ])

    plt.figure(figsize=(8, 6))
    sns.boxplot(x='Category', y='length', data=plot_data, palette='pastel')
    plt.title(f'Sequence Length Comparison for {model_name}', fontsize=16)
    plt.ylabel('Full Sequence Length')
    plt.xlabel('')
    plt.tight_layout()
    filename = os.path.join(output_dir, f'{model_name}_TP_vs_FN_length_comparison.png')
    plt.savefig(filename)
    plt.close()
    print(f"Saved length comparison plot to '{filename}'")

def get_aa_composition(sequences):
    """Calculates the average amino acid composition for a list of sequences."""
    compositions = []
    for seq in sequences:
        n_term = seq[:N_TERMINUS_LEN]
        if not n_term: continue # Skip empty sequences
        counts = Counter(n_term)
        total = len(n_term)
        compositions.append([counts.get(aa, 0) / total for aa in AMINO_ACIDS])
    if not compositions: # Handle case where all sequences were empty
        return pd.Series(0, index=AMINO_ACIDS)
    return pd.DataFrame(compositions, columns=AMINO_ACIDS).mean()

def plot_aa_composition_comparison(tp_df, fn_df, model_name, output_dir):
    """Compares and plots the N-terminal AA composition of TPs vs. FNs."""
    if tp_df.empty or fn_df.empty:
        return

    tp_comp = get_aa_composition(tp_df['Sequence'])
    fn_comp = get_aa_composition(fn_df['Sequence'])

    comp_df = pd.DataFrame({'True Positive': tp_comp, 'False Negative': fn_comp})

    comp_df.plot(kind='bar', figsize=(18, 7), width=0.8)
    plt.title(f'N-Terminus Amino Acid Composition: TP vs. FN ({model_name})', fontsize=16)
    plt.ylabel('Average Frequency')
    plt.xlabel('Amino Acid')
    plt.xticks(rotation=0)
    plt.legend()
    plt.tight_layout()
    filename = os.path.join(output_dir, f'{model_name}_TP_vs_FN_aa_composition.png')
    plt.savefig(filename)
    plt.close()
    print(f"Saved AA composition plot to '{filename}'")

def plot_physical_property_comparison(tp_df, fn_df, model_name, output_dir):
    """Compares N-terminal AA physical properties of TPs vs. FNs."""
    if tp_df.empty or fn_df.empty:
        return

    tp_comp = get_aa_composition(tp_df['Sequence'])
    fn_comp = get_aa_composition(fn_df['Sequence'])

    def sum_properties(composition):
        prop_sums = {}
        for prop, aas in AA_PHYSICAL_PROPERTIES.items():
            prop_sums[prop] = composition[aas].sum()
        return pd.Series(prop_sums)

    tp_props = sum_properties(tp_comp)
    fn_props = sum_properties(fn_comp)

    props_df = pd.DataFrame({'True Positive': tp_props, 'False Negative': fn_props})

    props_df.plot(kind='bar', figsize=(10, 6), width=0.8)
    plt.title(f'N-Terminus Physicochemical Properties: TP vs. FN ({model_name})', fontsize=16)
    plt.ylabel('Total Average Frequency')
    plt.xlabel('Property Group')
    plt.xticks(rotation=0)
    plt.legend()
    plt.tight_layout()
    filename = os.path.join(output_dir, f'{model_name}_TP_vs_FN_physical_properties.png')
    plt.savefig(filename)
    plt.close()
    print(f"Saved physical properties plot to '{filename}'")

# --- Main Execution ---

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
        return

    X_test_list = test_df['Sequence'].apply(create_augmented_features).tolist()
    X_test = np.array(X_test_list)
    X_test_scaled = scaler.transform(X_test)

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

    # --- 7. Visual Error Analysis for BOTH Models ---
    print("\n[6] Generating plots for visual error analysis...")
    # Create a directory to save plots
    output_dir = "error_analysis_plots"
    os.makedirs(output_dir, exist_ok=True)

    # Loop through both models to generate a full set of plots for each
    for model_to_analyze in ['svm', 'von_heijne']:
        print(f"\n--- Generating visualizations for the {model_to_analyze.upper()} model ---")

        pred_col = f'predicted_{model_to_analyze}'
        is_fp = (test_df_extended['label'] == 0) & (test_df_extended[pred_col] == 1)
        is_fn = (test_df_extended['label'] == 1) & (test_df_extended[pred_col] == 0)
        is_tp = (test_df_extended['label'] == 1) & (test_df_extended[pred_col] == 1)

        fp_df = test_df_extended[is_fp].copy()
        fn_df = test_df_extended[is_fn].copy()
        tp_df = test_df_extended[is_tp].copy()

        # Generate all plots for the current model
        plot_taxonomy_pie_charts(fp_df, fn_df, model_to_analyze.upper(), output_dir)
        plot_sequence_logos(tp_df, fn_df, model_to_analyze.upper(), output_dir)
        plot_length_comparison(tp_df, fn_df, model_to_analyze.upper(), output_dir)
        plot_aa_composition_comparison(tp_df, fn_df, model_to_analyze.upper(), output_dir)
        plot_physical_property_comparison(tp_df, fn_df, model_to_analyze.upper(), output_dir)

    print(f"\nAll plots have been saved to the '{output_dir}/' directory.")
    print("\n--- Analysis Complete ---")

if __name__ == '__main__':
    main()
