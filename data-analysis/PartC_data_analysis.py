import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter
import os


def compute_aa_freq(seq, aa_keys):
    """
    Compute normalized amino acid frequency (%) for a given sequence.
    Returns a dictionary {amino_acid: frequency}.
    """
    seq = seq.upper()
    counts = Counter(seq)
    total = sum(counts.values())
    if total == 0:
        return {aa: 0 for aa in aa_keys}
    return {aa: (counts.get(aa, 0) / total) * 100 for aa in aa_keys}


def average_freq(df_subset, aa_keys):
    """
    Calculate the average amino acid frequency across a dataframe subset.
    """
    freq_dicts = df_subset['aa_freq'].tolist()
    if not freq_dicts:
        return {aa: 0 for aa in aa_keys}
    avg = {aa: sum(d[aa] for d in freq_dicts) / len(freq_dicts) for aa in aa_keys}
    return avg


def main(input_data = "merged_dataset_with_seqs.tsv", output_path="analysis_output"):
    """
    Main function to run the data analysis and visualization script.
    """
    # --- Configuration ---
    # General output path for all generated files
    output_path = "analysis_output"
    
    # --- Setup ---
    # Define the path for plots within the general output directory
    plots_path = os.path.join(output_path, 'plots')
    
    # Create the output directories if they don't exist
    os.makedirs(plots_path, exist_ok=True)
    print(f"Output will be saved to the '{output_path}' directory.")

    # Set the global theme for all plots
    sns.set_theme(
        context="notebook",
        style="darkgrid",
        palette="muted",
        font="sans-serif",
        font_scale=1,
        color_codes=True,
        rc=None
    )

    # Load the dataset
    try:
        df = pd.read_csv(input_data, sep="\t")
    except FileNotFoundError:
        print("Error: The file '/content/merged_dataset_with_seqs.tsv' was not found.")
        print("Please make sure the dataset file is in the correct directory.")
        return

    # Create data subsets. Using .copy() to avoid SettingWithCopyWarning later.
    df_train = df[df["Type"] == "train"].copy()
    df_test = df[df["Type"] == "test"].copy()
    df_pos = df[df["label"] == 1].copy()
    df_neg = df[df["label"] == 0].copy()

    print("Data loaded successfully. Here are the first 5 rows:")
    print(df.head())

    # --- Analysis of Protein Lengths ---

    ## Plot 1: Density of protein lengths (Positive vs. Negative) in the whole dataset
    sns.set_style("whitegrid")
    plt.figure(figsize=(10, 6))
    sns.kdeplot(data=df, x='Protein_Length', hue='label', fill=True, common_norm=False, palette={1: 'lightgreen', 0: 'orangered'})
    plt.title('Density Plot of Protein Lengths In The Whole Dataset', fontsize=16)
    plt.xlabel('Protein Length (amino acids)', fontsize=12)
    plt.ylabel('Density', fontsize=12)
    plt.xlim(0, 3000)
    plt.tight_layout()
    plt.savefig(os.path.join(plots_path, 'Density_Plot_of_Protein_Lengths_Whole_Dataset.png'))
    

    ## Plot 2: Density of protein lengths in Train vs. Test sets
    plt.figure(figsize=(10, 6))
    sns.kdeplot(data=df, x='Protein_Length', hue='Type', fill=True, common_norm=False, palette={"train": 'lightgreen', "test": 'orangered'})
    plt.title('Density Plot of Protein Lengths in Train and Test sets', fontsize=16)
    plt.xlabel('Protein Length (amino acids)', fontsize=12)
    plt.ylabel('Density', fontsize=12)
    plt.xlim(0, 3000)
    plt.tight_layout()
    plt.savefig(os.path.join(plots_path, 'Density_Plot_of_Protein_Lengths_Train_vs_Test.png'))
    

    ## Plot 3: Density of protein lengths in the Train set (Positive vs. Negative)
    plt.figure(figsize=(10, 6))
    sns.kdeplot(data=df_train, x='Protein_Length', hue='label', fill=True, common_norm=False, palette={1: 'lightgreen', 0: 'orangered'})
    plt.title('Density Plot of Protein Lengths in Train Set', fontsize=16)
    plt.xlabel('Protein Length (amino acids)', fontsize=12)
    plt.ylabel('Density', fontsize=12)
    plt.xlim(0, 3000)
    plt.tight_layout()
    plt.savefig(os.path.join(plots_path, 'Density_Plot_of_Protein_Lengths_in_Train_Set.png'))
    

    ## Plot 4: Density of protein lengths in the Test set (Positive vs. Negative)
    plt.figure(figsize=(10, 6))
    sns.kdeplot(data=df_test, x='Protein_Length', hue='label', fill=True, common_norm=False, palette={1: 'lightgreen', 0: 'orangered'})
    plt.title('Density Plot of Protein Lengths in Test Set', fontsize=16)
    plt.xlabel('Protein Length (amino acids)', fontsize=12)
    plt.ylabel('Density', fontsize=12)
    plt.xlim(0, 3000)
    plt.tight_layout()
    plt.savefig(os.path.join(plots_path, 'Density_Plot_of_Protein_Lengths_in_Test_Set.png'))
    

    ## Plot 5: Boxplot of protein lengths by label
    plt.figure(figsize=(10, 6))
    sns.boxplot(data=df, x="label", y="Protein_Length", hue="label", palette={1: 'springgreen', 0: 'orangered'})
    plt.title('Boxplot of Protein Lengths by Label', fontsize=16)
    plt.xlabel('Label', fontsize=12)
    plt.ylabel('Protein Length (amino acids)', fontsize=12)
    plt.legend(title='Type', labels=['Negative', 'Positive'])
    plt.ylim(0, 3000)
    plt.tight_layout()
    plt.savefig(os.path.join(plots_path, 'Boxplot_of_Protein_Lengths_by_Train_Test.png'))
    

    # --- Analysis of Signal Peptide (SP) Lengths ---

    ## Plot 6: Distribution of SP cleavage sites
    plt.figure(figsize=(10, 6))
    df_positive_only = df[df['label'] == 1]
    sns.histplot(data=df_positive_only, x='Cleavage_Site', hue="Type", kde=True)
    plt.title('Distribution of SP Cleavage Sites (Train vs Test)', fontsize=16)
    plt.xlabel('SP Cleavage Site Position', fontsize=12)
    plt.ylabel('Count', fontsize=12)
    plt.tight_layout()
    plt.savefig(os.path.join(plots_path, 'Distribution_of_SP_Cleavage_Sites.png'))
    

    # --- Amino Acid Composition Analysis ---

    # Background amino acid frequencies from SwissProt/ExPASy
    aa_freq_expasy = {
        'A': 8.25, 'Q': 3.93, 'L': 9.64, 'S': 6.65, 'R': 5.52, 'E': 6.71,
        'K': 5.79, 'T': 5.36, 'N': 4.06, 'G': 7.07, 'M': 2.41, 'W': 1.1,
        'D': 5.46, 'H': 2.27, 'F': 3.86, 'Y': 2.92, 'C': 1.38, 'I': 5.9,
        'P': 4.74, 'V': 6.85
    }
    aas = list(aa_freq_expasy.keys())

    ## Plot 7: Comparative AA composition (Positive vs. Negative)
    df_pos['aa_freq'] = df_pos['Sequence'].apply(lambda seq: compute_aa_freq(seq, aas))
    df_neg['aa_freq'] = df_neg['Sequence'].apply(lambda seq: compute_aa_freq(seq, aas))
    
    avg_pos = average_freq(df_pos, aas)
    avg_neg = average_freq(df_neg, aas)

    pos_freqs = [avg_pos[aa] for aa in aas]
    neg_freqs = [avg_neg[aa] for aa in aas]
    bg_freqs  = [aa_freq_expasy[aa] for aa in aas]

    x = np.arange(len(aas))
    width = 0.35
    plt.figure(figsize=(14, 6))
    plt.bar(x - width/2, pos_freqs, width, label='Positive (label=1)', alpha=0.7, color='lightgreen')
    plt.bar(x + width/2, neg_freqs, width, label='Negative (label=0)', alpha=0.7, color='salmon')
    plt.plot(x, bg_freqs, color='black', marker='o', linestyle='--', label='Background (ExPASy)')
    plt.xticks(x, aas, fontsize=10)
    plt.xlabel("Amino Acid")
    plt.ylabel("Frequency (%)")
    plt.title("Comparative Amino Acid Composition: Positive vs Negative")
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(plots_path, 'AA_Composition_Positive_vs_Negative.png'))
    

    ## Plot 8: Comparative AA composition (Train vs. Test)
    df_train['aa_freq'] = df_train['Sequence'].apply(lambda seq: compute_aa_freq(seq, aas))
    df_test['aa_freq'] = df_test['Sequence'].apply(lambda seq: compute_aa_freq(seq, aas))

    avg_train = average_freq(df_train, aas)
    avg_test = average_freq(df_test, aas)

    train_freqs = [avg_train[aa] for aa in aas]
    test_freqs = [avg_test[aa] for aa in aas]

    plt.figure(figsize=(14, 6))
    plt.bar(x - width/2, train_freqs, width, label='Train', alpha=0.7, color='plum')
    plt.bar(x + width/2, test_freqs, width, label='Test', alpha=0.7, color='skyblue')
    plt.plot(x, bg_freqs, color='black', marker='o', linestyle='--', label='Background (ExPASy)')
    plt.xticks(x, aas, fontsize=10)
    plt.xlabel("Amino Acid", fontsize=12)
    plt.ylabel("Frequency (%)", fontsize=12)
    plt.title("Comparative Amino Acid Composition: Train vs. Test", fontsize=14, fontweight='bold')
    plt.legend(loc='upper right')
    plt.tight_layout()
    plt.savefig(os.path.join(plots_path, 'AA_Composition_Train_vs_Test.png'))
    

    # --- Taxonomic Classification Analysis ---

    ## Plot 9: Bar chart of Kingdom distribution
    plt.figure(figsize=(10, 6))
    sns.countplot(data=df, x='Kingdom', hue='Type', palette={"train": 'plum', "test": 'skyblue'})
    plt.title('Taxonomy Classification - Train vs Test', fontsize=16)
    plt.xlabel('Kingdom', fontsize=12)
    plt.ylabel('Count', fontsize=12)
    plt.tight_layout()
    plt.savefig(os.path.join(plots_path, 'Taxonomy_Classification_Bar_Chart.png'))
    

    ## Plot 10: Pie charts of Kingdom distribution
    df_train_kingdom_n = df_train['Kingdom'].value_counts()
    df_test_kingdom_n = df_test['Kingdom'].value_counts()
    colors = ["salmon", "lightblue", "lightgreen", "hotpink"]

    fig, axes = plt.subplots(1, 2, figsize=(14, 7))
    axes[0].pie(df_train_kingdom_n, labels=df_train_kingdom_n.index, colors=colors, autopct='%1.1f%%', startangle=90)
    axes[0].set_title('Train Set: Distribution of Kingdoms')
    axes[0].axis('equal')
    axes[1].pie(df_test_kingdom_n, labels=df_test_kingdom_n.index, colors=colors, autopct='%1.1f%%', startangle=90)
    axes[1].set_title('Test Set: Distribution of Kingdoms')
    axes[1].axis('equal')
    plt.savefig(os.path.join(plots_path, 'Kingdom_Classification_Pie_Charts.png'))
    
    
    # Plot 11: Pie charts of Organism distribution
    df_train_processed = df_train.copy()
    df_test_processed = df_test.copy()

    # 1. Initial Organism Counts (using original df_train for the counts)
    # Note: df is not modified, so we can keep the first line as is for labels.
    labels = df['Organism'].unique()
    df_train_organism_n = df_train_processed['Organism'].value_counts()
    df_test_organism_n = df_test_processed['Organism'].value_counts()
    colors = ["salmon", "lightblue", "lightgreen", "hotpink"]

    # set a threshold and change the name of the organisms below that to "Other"
    threshold = 500 

    organisms_to_group = df_train_organism_n[df_train_organism_n < threshold].index

    # 2. Apply grouping to the COPIED training set
    df_train_processed.loc[df_train_processed['Organism'].isin(organisms_to_group), 'Organism'] = 'Other'


    organisms_to_keep_train = df_train_processed['Organism'].unique() # Get the new labels from the modified train set

    # 3. Apply grouping to the COPIED test set
    # Group anything in the test set that is NOT one of the 'kept' train organisms
    df_test_processed.loc[~df_test_processed['Organism'].isin(organisms_to_keep_train), 'Organism'] = 'Other'

    # 4. Re-calculate value counts using the NEW processed dataframes
    df_train_organism_n_new = df_train_processed['Organism'].value_counts()
    df_test_organism_n_new = df_test_processed['Organism'].value_counts()

    df_train_organism_n_new


    ### Plotting
    fig, axes = plt.subplots(1, 2, figsize=(14, 7))
    # Plot for positive data
    axes[0].pie(df_train_organism_n_new, labels=df_train_organism_n_new.index, colors=colors, autopct='%1.1f%%', startangle=90)
    axes[0].set_title('Train Set: Distribution of Organisms')
    axes[0].axis('equal')

    # Plot for negative data
    axes[1].pie(df_test_organism_n_new, labels=df_test_organism_n_new.index, colors=colors, autopct='%1.1f%%', startangle=90)
    axes[1].set_title('Test Set: Distribution of Organisms')
    axes[1].axis('equal')

    #save
    plt.savefig(os.path.join(plots_path, 'Organism_Classification_Pie_Charts.png'))
    

    # --- Sequence Motif Extraction ---

    motifs = []
    for _, row in df_pos.iterrows():
        sequence = row['Sequence']
        cleavage_site = row['Cleavage_Site']
        
        # Ensure cleavage site is a valid number and within sequence bounds
        try:
            # Convert cleavage site to 0-based index
            cleavage_site_idx = int(cleavage_site) - 1
            start = cleavage_site_idx - 13
            end = cleavage_site_idx + 2
            if start >= 0 and end <= len(sequence):
                motifs.append(sequence[start:end])
        except (ValueError, TypeError):
            # Skip if cleavage site is not a valid number
            continue

    # Define the full path for the motifs file
    motifs_filepath = os.path.join(output_path, 'motifs.txt')
    
    # Write the extracted motifs to the text file
    with open(motifs_filepath, 'w') as f:
        for motif in motifs:
            f.write(f"{motif}\n")
            
    print(f"\nSuccessfully extracted and saved {len(motifs)} motifs to '{motifs_filepath}'")


if __name__ == '__main__':
    main()