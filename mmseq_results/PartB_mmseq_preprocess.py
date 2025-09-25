import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split, StratifiedKFold
import subprocess
import os
import tabulate

def parse_fasta(file_path):
    """
    Parses a FASTA file and returns a dictionary of sequences.
    """
    sequences = {}
    current_id = None
    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                current_id = line[1:].split()[0]
                sequences[current_id] = ''
            elif current_id:
                sequences[current_id] += line
    return sequences

def main():
    """
    Main function to execute the data processing and analysis pipeline.
    """
    # Define file paths
    # Please update these paths to the correct location of your files
    neg_cluster_file = "neg_cluster-results_cluster.tsv"
    pos_cluster_file = "pos_cluster-results_cluster.tsv"
    og_neg_file = "negative_dataset.tsv"
    og_pos_file = "positive_dataset_sp_cleavage.tsv"
    neg_fasta_file = "negative_dataset.fasta"
    pos_fasta_file = "positive_dataset_sp_cleavage.fasta"

    # Check if input files exist
    for f in [neg_cluster_file, pos_cluster_file, og_neg_file, og_pos_file, neg_fasta_file, pos_fasta_file]:
        if not os.path.exists(f):
            print(f"Error: Input file not found at {f}")
            return

    # Load cluster data
    neg_cluster = pd.read_csv(neg_cluster_file, sep="\t", header=None)
    pos_cluster = pd.read_csv(pos_cluster_file, sep="\t", header=None)

    # Load original datasets
    og_neg = pd.read_csv(og_neg_file, sep="\t")
    og_pos = pd.read_csv(og_pos_file, sep="\t")

    # Get representative sequences from clusters
    rep_neg = neg_cluster[0].unique()
    rep_pos = pos_cluster[0].unique()

    # Filter original datasets to keep only representative sequences
    og_neg_filter = og_neg[og_neg['Accession'].isin(rep_neg)].reset_index(drop=True)
    og_pos_filter = og_pos[og_pos['Accession'].isin(rep_pos)].reset_index(drop=True)

    # Save filtered datasets
    og_neg_filter.to_csv("filtered_neg_dataset.tsv", sep="\t", index=False)
    og_pos_filter.to_csv("filtered_pos_dataset.tsv", sep="\t", index=False)

    # Create training and testing splits
    pos_train, pos_test = train_test_split(og_pos_filter, test_size=0.2, random_state=42)
    neg_train, neg_test = train_test_split(og_neg_filter, test_size=0.2, random_state=42)

    # Add labels
    pos_train['label'] = 1
    neg_train['label'] = 0
    pos_test['label'] = 1
    neg_test['label'] = 0

    # Combine and shuffle datasets
    train_df = pd.concat([pos_train, neg_train], axis=0).sample(frac=1, random_state=42).reset_index(drop=True)
    test_df = pd.concat([pos_test, neg_test], axis=0).sample(frac=1, random_state=42).reset_index(drop=True)

    # Create 5-fold cross-validation subsets
    train_df['fold'] = -1
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    for fold_num, (train_index, val_index) in enumerate(skf.split(train_df, train_df['label'])):
        train_df.loc[val_index, 'fold'] = fold_num

    print("Value counts for each fold:")
    print(train_df['fold'].value_counts())
    print("\nVerifying the positive/negative ratio in each fold:")
    print(train_df.groupby('fold')['label'].value_counts(normalize=True))

    # Add 'Type' column and merge train and test sets
    train_df["Type"] = "train"
    test_df["Type"] = "test"
    merged_df = pd.concat([train_df, test_df], axis=0)

    # Save merged dataset
    merged_df.to_csv("merged_dataset.tsv", sep="\t", index=False)

    # Merge FASTA files
    merged_fasta_file = "merged.fasta"
    with open(merged_fasta_file, 'wb') as outfile:
        for fasta_file in [neg_fasta_file, pos_fasta_file]:
            with open(fasta_file, 'rb') as infile:
                outfile.write(infile.read())

    # Parse FASTA and add sequences to the merged dataframe
    fasta_sequences = parse_fasta(merged_fasta_file)
    merged_df["Sequence"] = merged_df['Accession'].map(fasta_sequences)

    # Save the final merged dataset with sequences
    merged_df.to_csv("merged_dataset_with_seqs.tsv", sep="\t", index=False)

    # --- Generate Summary Table ---
    train_n = merged_df[merged_df["Type"] == "train"].shape[0]
    test_n = merged_df[merged_df["Type"] == "test"].shape[0]
    train_label_n_1 = merged_df[(merged_df["Type"] == "train") & (merged_df["label"] == 1)].shape[0]
    train_label_n_0 = merged_df[(merged_df["Type"] == "train") & (merged_df["label"] == 0)].shape[0]
    test_label_n_1 = merged_df[(merged_df["Type"] == "test") & (merged_df["label"] == 1)].shape[0]
    test_label_n_0 = merged_df[(merged_df["Type"] == "test") & (merged_df["label"] == 0)].shape[0]

    neg_before_cluster_n = og_neg.shape[0]
    pos_before_cluster_n = og_pos.shape[0]
    neg_after_cluster_n = len(rep_neg)
    pos_after_cluster_n = len(rep_pos)

    total_negative = train_label_n_0 + test_label_n_0
    total_positive = train_label_n_1 + test_label_n_1

    summary_data = [
        {"Dataset": "Train", "Negative": train_label_n_0, "Positive": train_label_n_1},
        {"Dataset": "Test", "Negative": test_label_n_0, "Positive": test_label_n_1},
        {"Dataset": "Total", "Negative": total_negative, "Positive": total_positive},
    ]

    summary_df = pd.DataFrame(summary_data)

    print("\n--- Data Summary ---")
    print(summary_df.to_markdown(index=False))

    print(f"\nNegative samples before clustering: {neg_before_cluster_n}")
    print(f"Positive samples before clustering: {pos_before_cluster_n}")
    print(f"Negative samples after clustering: {neg_after_cluster_n}")
    print(f"Positive samples after clustering: {pos_after_cluster_n}")


if __name__ == "__main__":
    main()