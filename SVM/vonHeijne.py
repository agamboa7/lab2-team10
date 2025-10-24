from typing import final
import pandas as pd
import numpy as np
from sklearn.metrics import precision_recall_curve, f1_score, precision_score, recall_score, matthews_corrcoef
import matplotlib.pyplot as plt
from numpy import trapezoid
import seaborn as sns
from tabulate import tabulate



# Constants
WINDOW_START = -13
WINDOW_END = 2
WINDOW_SIZE = WINDOW_END - WINDOW_START + 1
PSEUDOCOUNT = 1
SCAN_REGION = 90
AMINO_ACIDS = 'ACDEFGHIKLMNPQRSTVWY'
aa_freq_expasy = {
    'A': 8.25, 'R': 5.52, 'N': 4.06, 'D': 5.46, 'C': 1.38, 'Q': 3.93,
    'E': 6.71, 'G': 7.07, 'H': 2.27, 'I': 5.9,  'L': 9.64, 'K': 5.79,
    'M': 2.41, 'F': 3.86, 'P': 4.74, 'S': 6.65, 'T': 5.36, 'W': 1.1,
    'Y': 2.92, 'V': 6.85
}
background_frequencies = {aa: freq / 100.0 for aa, freq in aa_freq_expasy.items()}


# Rest of the thing 💠 🏮 🍨 🛡 📞 💾 📈 💐 🍉 ▶️ 🛢 🐢 ⚰ 👼 🎷 🚾 🕑 👪 🌬 🔊 😄 👅 🎭 🚟 🦃 🐃 🐳 🍩 🎰 ✌️

def predict_scores(data_df, score_matrix):
    """Calculates scores for each sequence."""
    all_scores = []
    for _, row in data_df.iterrows():
        seq = row["Sequence"]
        max_score = float('-inf')
        scan_limit = min(len(seq), SCAN_REGION) - WINDOW_SIZE + 1

        for i in range(scan_limit):
            window = seq[i : i + WINDOW_SIZE]
            if len(window) == WINDOW_SIZE and all(aa in AMINO_ACIDS for aa in window):
                current_score = sum(score_matrix[aa][pos] for pos, aa in enumerate(window))
                if current_score > max_score:
                    max_score = current_score
        all_scores.append(max_score if max_score != float('-inf') else -999)
    return all_scores


def plot_pswm_heatmap(pswm):
    """
    Generates a heatmap of the PSWM with high (positive) values as green
    and low (negative) values as red.
    """
    print("\n5. Generating and saving PSWM Heatmap...")

    pswm_df = pd.DataFrame(pswm).T
    window_cols = [f'P{i:+d}' for i in range(WINDOW_START, WINDOW_END + 1)]
    pswm_df.columns = window_cols

    vmax = pswm_df.abs().max().max()
    vmin = -vmax

    plt.figure(figsize=(WINDOW_SIZE * 0.5, len(AMINO_ACIDS) * 0.5)) 

    # 3. Create the heatmap
    sns.heatmap(
        pswm_df,
        cmap='RdYlGn', 
        vmin=vmin,
        vmax=vmax,
        center=0, 
        annot=True,
        fmt=".2f",
        linewidths=.5,
        cbar_kws={'label': 'Log-Odds Score (PSWM)'}
    )

    plt.title('Position-Specific Weight Matrix (PSWM) Heatmap', fontsize=14)
    plt.xlabel('Position Relative to Cleavage Site', fontsize=12)
    plt.ylabel('Amino Acid', fontsize=12)

    plt.tight_layout() 
    plt.savefig('pswm_heatmap.png', dpi=300)
    print("Heatmap saved as 'pswm_heatmap.png'")
    plt.show()



def training(train_df, val_df):
    """
    Trains the PSWM and finds the best threshold.

    """
    # Get all the motifs from positive seqs to train matrix
    motifs = []
    for _, row in train_df[train_df["label"] == 1].iterrows():
        cs = int(row['Cleavage_Site']) - 1
        motif = row["Sequence"][cs + WINDOW_START : cs + WINDOW_END + 1]
        if len(motif) == WINDOW_SIZE:
            motifs.append(motif)

    #Create Matrix
    score_matrix = {aa: [PSEUDOCOUNT] * WINDOW_SIZE for aa in AMINO_ACIDS}
    for motif in motifs:
        for i, aa in enumerate(motif):
            if aa in score_matrix: score_matrix[aa][i] += 1

    n = len(motifs)
    divisor = n + len(AMINO_ACIDS)
    pspm = {aa: [count / divisor for count in score_matrix[aa]] for aa in AMINO_ACIDS}
    pswm = {aa: [np.log(pspm[aa][i] / background_frequencies[aa]) for i in range(WINDOW_SIZE)] for aa in AMINO_ACIDS}

    ### Finding a threshold
    y_validation_true = val_df["label"].tolist()
    y_validation_scores = predict_scores(val_df, pswm)

    precision, recall, thresholds = precision_recall_curve(y_validation_true, y_validation_scores)
    fscore = (2 * precision * recall) / (precision + recall + 1e-9)
    optimal_threshold = thresholds[np.argmax(fscore)]

    #Return the precision, recall for later plotting
    return pswm, optimal_threshold, precision, recall


def prediction(test_df, score_matrix, threshold):
    """Predicts labels for the test set."""
    test_scores = predict_scores(test_df, score_matrix)
    return [1 if score >= threshold else 0 for score in test_scores]


def plot_averaged_pr_curve(all_precisions, all_recalls):
    """
    Takes lists of precision and recall arrays from a CV run and plots the averaged curve.
    """
    print("\n3. Generating and saving averaged Precision-Recall curve plot...")

    # Interpolate all curves to a common recall axis
    mean_recall_axis = np.linspace(0, 1, 100)
    interpolated_precisions = []
    for recall, precision in zip(all_recalls, all_precisions):
        recall_increasing = recall[::-1]
        precision_increasing = precision[::-1]
        interp_precision = np.interp(mean_recall_axis, recall_increasing, precision_increasing)
        interpolated_precisions.append(interp_precision)

    interpolated_precisions = np.array(interpolated_precisions)

    # Calculate mean and standard deviation
    mean_precision_curve = np.mean(interpolated_precisions, axis=0)
    std_precision = np.std(interpolated_precisions, axis=0)
    pr_auc = trapezoid(mean_precision_curve, mean_recall_axis)

    # Generate the plot using your specified design
    plt.figure(figsize=(8, 6))
    plt.plot(mean_recall_axis, mean_precision_curve, color='b',
              label=fr'Mean PR Curve (AUC = {pr_auc:.2f})')

    tprs_upper = np.minimum(mean_precision_curve + std_precision, 1)
    tprs_lower = np.maximum(mean_precision_curve - std_precision, 0)
    plt.fill_between(mean_recall_axis, tprs_lower, tprs_upper, color='grey', alpha=.2,
                      label=r'$\pm$ 1 std. dev.')

    plt.xlim([-0.05, 1.05])
    plt.ylim([-0.05, 1.05])
    plt.xlabel('Recall', fontsize=12)
    plt.ylabel('Precision', fontsize=12)
    plt.title('Averaged Precision-Recall Curve from Cross-Validation', fontsize=14)
    plt.legend(loc="lower left")
    plt.grid(True)
    plt.savefig('averaged_precision_recall_curve.png', dpi=300)
    print("Plot saved as 'averaged_precision_recall_curve.png'")
    plt.show()

def main(df):
    """
    Runs cross-validation, reports performance, and generates the averaged PR plot.
    """
    all_f1_scores, all_precision_scores, all_recall_scores, all_mcc_scores = [], [], [], []

    # Vars for plot
    pr_curve_precisions = []
    pr_curve_recalls = []
    all_thres = []

    print("1. Starting 5-fold cross-validation...")
    for j in range(5):
        test_index, validation_index = j, (j + 1) % 5
        training_indices = [k for k in range(5) if k not in {test_index, validation_index}]

        print(f"  - Run {j+1}/5: Test Fold: {test_index}, Validation Fold: {validation_index}, Training Folds: {training_indices}")

        test_df = df[df["fold"] == test_index]
        validation_df = df[df["fold"] == validation_index]
        training_df = df[df["fold"].isin(training_indices)]

        pswm, threshold, fold_precision, fold_recall = training(training_df, validation_df)
        all_thres.append(threshold)

        pr_curve_precisions.append(fold_precision)
        pr_curve_recalls.append(fold_recall)

        y_pred = prediction(test_df, pswm, threshold)
        y_true = test_df["label"].tolist()

        all_f1_scores.append(f1_score(y_true, y_pred, zero_division=0))
        all_precision_scores.append(precision_score(y_true, y_pred, zero_division=0))
        all_recall_scores.append(recall_score(y_true, y_pred, zero_division=0))
        all_mcc_scores.append(matthews_corrcoef(y_true, y_pred))

    ### Calculate the metrics
    #f1
    mean_f1 = np.mean(all_f1_scores)
    se_f1 = np.std(all_f1_scores) / np.sqrt(5)
    final_f1 = f"{mean_f1:.4f} ± {se_f1:.4f}"

    #precision
    mean_precision = np.mean(all_precision_scores)
    se_precision = np.std(all_precision_scores) / np.sqrt(5)
    final_precision = f"{mean_precision:.4f} ± {se_precision:.4f}"

    #recall
    mean_recall = np.mean(all_recall_scores)
    se_recall = np.std(all_recall_scores) / np.sqrt(5)
    final_recall = f"{mean_recall:.4f} ± {se_recall:.4f}"

    #mcc
    mean_mcc = np.mean(all_mcc_scores)
    se_mcc = np.std(all_mcc_scores) / np.sqrt(5)
    final_mcc = f"{mean_mcc:.4f} ± {se_mcc:.4f}"

    #avg threshold
    mean_threshold = np.mean(all_thres)
    se_threshold = np.std(all_thres) / np.sqrt(5)
    final_threshold = f"{mean_threshold:.4f} ± {se_threshold:.4f}"

    table_data = [
        ("F1 Score (Mean ± SE)", final_f1),
        ("Precision (Mean ± SE)", final_precision),
        ("Recall (Mean ± SE)", final_recall),
        ("Avg Threshold (Mean ± SE)", final_threshold),
        ("MCC (Mean ± SE)", final_mcc),
    ]
    
    # 2. Define the headers
    headers = ["Metric", "Value"]

    # 3. Print the title and the table
    print("\n2. Cross-Validation Final Performance Metrics:")

    # 'grid' format provides a clean, box-style table
    print(tabulate(table_data, headers=headers, tablefmt="github"))
    
    #save the table
    with open("cross_validation_results.txt", "w") as f:
        f.write(tabulate(table_data, headers=headers, tablefmt="github"))

    # PLot PR curve
    plot_averaged_pr_curve(pr_curve_precisions, pr_curve_recalls)

    #Plot the heatmap of matrix
    plot_pswm_heatmap(pswm)

if __name__ == "__main__":
  input_tsv = "merged_dataset_with_seqs.tsv"

  df = pd.read_csv(input_tsv, sep="\t")
  main(df)