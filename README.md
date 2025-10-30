# lab2-team10
This repository contains material relevant to the Laboratory of Bioinformatics 2 course, which is part of the Bioinformatics Master's Degree at the University of Bologna. 

Team Members:
Andrea Arriola Gamboa, Kristina Djikic, Başak Akkoyun, Betül Yalçın, Deniz Ertuğrul



# Signal Peptide Detection: Machine Learning Pipeline

This repository contains a comprehensive machine learning pipeline for predicting signal peptides (SP) in eukaryotic protein sequences. The project compares multiple approaches including the classic von Heijne algorithm, Support Vector Machines (SVM), and Deep Learning models.

---

## Table of Contents

1. [Project Overview](#project-overview)
2. [Repository Structure](#repository-structure)
3. [Pipeline Architecture](#pipeline-architecture)
4. [Key Results](#key-results)
5. [Detailed Module Descriptions](#detailed-module-descriptions)
6. [Running the Pipeline](#running-the-pipeline)
7. [Dependencies](#dependencies)

---

## Project Overview

### What are Signal Peptides?

Signal peptides are short amino acid sequences (typically 14-30 residues) located at the N-terminus of many eukaryotic proteins. They direct nascent proteins to the endoplasmic reticulum for secretion or membrane insertion, and are cleaved off during translocation. Accurately predicting signal peptides is crucial for:

- Understanding protein localization
- Characterizing secretome composition
- Identifying therapeutic targets
- Distinguishing transmembrane proteins from secreted proteins

### Project Goals

1. Develop a comprehensive dataset of positive (proteins with experimentally verified SPs) and negative (proteins without SPs) examples
2. Apply multiple machine learning approaches to SP classification
3. Compare performance, particularly on challenging cases like transmembrane proteins
4. Create generalizable models for broad eukaryotic protein families

---

## Repository Structure

```
lab2-team10/
├── README.md
│
├── data-collection/          # Phase A: Data acquisition
│   ├── dataset_gathering.py  # UniProt API data collection
│   ├── positive_dataset_sp_cleavage.tsv
│   ├── negative_dataset.tsv
│   └── README.md
│
├── mmseq_results/            # Phase B: Clustering & preprocessing
│   ├── PartB_mmseq_preprocess.py
│   ├── merged_dataset.tsv
│   ├── merged_dataset_with_seqs.tsv
│   └── README.md
│
├── data-analysis/            # Phase C: Exploratory data analysis
│   ├── PartC_data_analysis.py
│   ├── plots/
│   └── README.md
│
├── model-creation/           # Phase D: Von Heijne algorithm
│   ├── model_creation.py
│   ├── cross_validation_results.txt
│   └── README.md
│
├── SVM/                      # Phase E: SVM with feature selection
│   ├── SVM.py
│   ├── vonHeijne.py (for comparison)
│   ├── model_comparison.py
│   ├── plots/
│   └── README.md
│
└── deep_learning/           # Phase F: Deep learning (CNN-BiLSTM)
    ├── deep_learning.py
    ├── complete_comparison.py
    ├── vonHeijne.py (for comparison)
    ├── SVM.py (for comparison)
    └── README.md
```

---

## Pipeline Architecture

### End-to-End Workflow

```
Phase A: Data Collection
    ↓
    └─→ Fetch positive/negative datasets from UniProt
        • Positive: Proteins with experimental SP evidence
        • Negative: Secreted proteins without SPs
        • Output: 2,932 positive + 20,615 negative sequences

Phase B: Data Preprocessing & Clustering
    ↓
    └─→ Remove sequence redundancy via MMseqs2 clustering
        • 90% sequence identity threshold
        • Select representative sequences
        • Output: 1,093 positive + 8,934 negative (non-redundant)

Phase C: Exploratory Data Analysis
    ↓
    └─→ Characterize datasets
        • Protein length distributions
        • SP length distributions
        • Amino acid composition analysis
        • Taxonomic distributions
        • Cleavage site motif analysis

Phase D: Von Heijne Algorithm
    ↓
    └─→ Train Position-Specific Weight Matrix (PSWM)
        • Uses window [-13, +2] relative to cleavage site
        • Cross-validated on 5 folds
        • F1 Score: 0.724 ± 0.016

Phase E: SVM with Feature Engineering
    ↓
    └─→ Machine learning with advanced features
        • 20 amino acid composition features (N-terminus)
        • Physicochemical properties (hydrophobicity, charge, etc.)
        • Feature selection via RandomForest
        • GridSearchCV hyperparameter optimization
        • F1 Score: 0.86, MCC: 0.847

Phase F: Deep Learning (CNN + BiLSTM)
    ↓
    └─→ Deep neural network classifier
        • 1D CNN (kernel size 17) for local motif extraction
        • 2-layer Bidirectional LSTM for long-range dependencies
        • Early stopping on validation MCC
        • F1 Score: 0.93, MCC: 0.92
```

---


## Key Results

### Model Comparison Summary

| Metric | Von Heijne | SVM | Deep Learning |
|--------|-----------|-----|---------------|
| **Overall Accuracy** | 0.94 | 0.97 | **0.99** |
| **F1 Score (Positive)** | 0.71 | 0.86 | **0.93** |
| **Precision (Positive)** | 0.72 | 0.87 | **0.94** |
| **Recall (Positive)** | 0.70 | 0.85 | **0.93** |
| **MCC** | 0.694 | 0.847 | **0.92** |
| **FPR on TM Proteins** | 0.206 | 0.127 | **0.036** |

### Key Findings

1. **Deep Learning dominates:** CNN-BiLSTM achieves highest accuracy (99%) and best handling of transmembrane proteins (FPR: 3.6%)

2. **SVM strong middle-ground:** Good performance (97% accuracy) with fast training and interpretable features

3. **Von Heijne still competitive:** Traditional method achieves 94% accuracy, useful as baseline/reference

4. **TM protein challenge:** All models struggle more with TM proteins (higher false positive rate) due to similar N-terminal hydrophobicity. Deep learning best mitigates this.

5. **Feature importance:** Hydrophobicity and TM propensity are top predictive features, confirming biological understanding of SP structure

---

## Detailed Module Descriptions

### 1. **data-collection/** - Phase A: Data Gathering

**Purpose:** Fetch and preprocess raw protein data from UniProtKB

**Script:** `dataset_gathering.py`

**Key Features:**
- Robust API client with retry logic for handling server errors
- Pagination to handle large search results (~3000 positive, ~20,600 negative sequences)
- Advanced filtering in addition to UniProt search:
  - Positive: Only cleaved SPs ≥14 residues
  - Negative: Secreted proteins without any SP evidence
- Outputs both TSV (metadata) and FASTA (sequences) formats

**UniProt Queries:**
- **Positive:** `(fragment:false) AND (taxonomy_id:2759) AND (length:[40 TO *]) AND (reviewed:true) AND (existence:1) AND (ft_signal_exp:*)`
- **Negative:** Secreted proteins with cytosolic/nuclear/organellar localization, excluding those with signal peptides

**Outputs:**
- `positive_dataset_sp_cleavage.tsv` / `.fasta`
- `negative_dataset.tsv` / `.fasta`

**Kingdom Distribution:**
| | Metazoa | Plants | Fungi | Other |
|---|---------|--------|-------|-------|
| Positive | 82.5% | 10.6% | 5.6% | 1.2% |
| Negative | 60.2% | 19.9% | 18.1% | 1.7% |

---

### 2. **mmseq_results/** - Phase B: Clustering & Preprocessing

**Purpose:** Remove sequence redundancy and create train/test splits with cross-validation folds

**Key Steps:**

1. **MMseqs2 Clustering** (run on VM):
   ```bash
   mmseqs easy-cluster input.fa cluster-results tmp --min-seq-id 0.3 -c 0.4 --cov-mode 0
   ```
   - 30% minimum sequence identity
   - Reduces 23,547 → 10,027 representative sequences

2. **Representative Sequence Selection:**
   - Extract cluster representatives (first member = highest similarity to all cluster members)
   - Filter original datasets to keep only representatives
   - Final: 1,093 positive + 8,934 negative

3. **Train/Test Split:**
   - 80% training / 20% test stratified split
   - Ensures class balance in both sets

4. **5-Fold Stratified Cross-Validation:**
   - Training set divided into 5 balanced folds
   - Used for hyperparameter tuning and model validation

5. **Sequence Integration:**
   - Merge original FASTA sequences with metadata
   - Output: `merged_dataset_with_seqs.tsv` (complete dataset with sequences and annotations)

**Output Datasets:**

| Dataset | Negative | Positive | Total |
|---------|----------|----------|-------|
| Train | 7,147 | 874 | 8,021 |
| Test | 1,787 | 219 | 2,006 |
| **Total** | **8,934** | **1,093** | **10,027** |

---

### 3. **data-analysis/** - Phase C: Exploratory Data Analysis

**Purpose:** Characterize training and test datasets to understand their properties

**Script:** `PartC_data_analysis.py`

**Analyses Performed:**

1. **Protein Length Distribution:**
   - Positive vs. negative sequences
   - Train vs. test sets
   - Density plots and box plots reveal:
     - Median positive length: ~300 aa
     - Median negative length: ~450 aa
     - Negative sequences are longer on average

2. **Signal Peptide Length Distribution:**
   - Cleavage site position analysis
   - Histogram shows most SPs are 14-30 residues (typical range)

3. **Amino Acid Composition:**
   - Compare SP sequences vs. SwissProt background
   - Key observation: SPs are enriched in:
     - Small hydrophobic residues (Ala, Val, Leu)
     - Deplete in charged residues (Arg, Asp, Glu)

4. **Taxonomic Distribution:**
   - Kingdom-level breakdown (Metazoa, Plants, Fungi, Other)
   - Consistent distribution in train/test sets
   - Organism-level pie charts

5. **Cleavage Site Motif Analysis:**
   - Extract [-13, +2] windows around cleavage sites
   - Generate sequence logos using WebLogo
   - Reveals conserved patterns for manual inspection

**Key Finding:** Train and test sets have consistent distributions across all metrics, validating dataset split quality.

---

### 4. **model-creation/** - Phase D: Von Heijne Algorithm

**Purpose:** Traditional statistical approach to SP detection based on Position-Specific Weight Matrices

**Script:** `model_creation.py`

**Algorithm Overview:**

1. **Extract Training Motifs:**
   - From positive training proteins, extract [-13, +2] windows at cleavage sites
   - Collect all motifs from folds designated as training set

2. **Build Position-Specific Probability Matrix (PSPM):**
   ```
   PSPM[aa][pos] = (count[aa] at pos + pseudocount) / (total + background)
   ```
   - Pseudocount = 1 for Laplace smoothing
   - Prevents zero probabilities

3. **Compute Position-Specific Weight Matrix (PSWM):**
   ```
   PSWM[aa][pos] = log(PSPM[aa][pos] / background_freq[aa])
   ```
   - Background frequencies from SwissProt composition
   - Log-odds scoring

4. **Prediction on Test/Validation:**
   - Scan protein N-terminus (first 90 residues)
   - Score all windows of length 16
   - Maximum score = protein's SP score

5. **Threshold Optimization:**
   - Use validation set to find optimal threshold
   - F1-score maximization via Precision-Recall curve
   - Threshold ≈ 6.73

**Cross-Validation Performance:**

| Metric | Value |
|--------|-------|
| F1 Score | 0.724 ± 0.016 |
| Precision | 0.750 ± 0.014 |
| Recall | 0.706 ± 0.035 |
| MCC | 0.694 ± 0.016 |

**Outputs:**
- `averaged_precision_recall_curve.png` - PR curve from all 5 folds
- `pswm_heatmap.png` - Visualization of weight matrix

---

### 5. **SVM/** - Phase E: Support Vector Machine with Feature Engineering

**Purpose:** Apply advanced machine learning with rich biochemical features

**Scripts:** 
- `SVM.py` - Main SVM pipeline with feature engineering
- `vonHeijne.py` - For comparison
- `model_comparison.py` - Side-by-side evaluation

**Feature Engineering (26 total features):**

**1. Amino Acid Composition (20 features):**
- Frequency of each of the 20 standard amino acids in N-terminal region (22 residues)

**2. Physicochemical Properties (6 features):**
- **n-region charge** (first 5 residues): Sum of charges from basic/acidic residues
- **h-region hydrophobicity** (residues 5-16): Average Kyte-Doolittle hydrophobicity
- **h-region TM propensity** (residues 5-16): Transmembrane helix propensity
- **Total N-terminal hydrophobicity**: Average KD across all 22 residues
- **Total N-terminal TM propensity**: Average TM across all 22 residues
- **Total N-terminal α-helix propensity**: Average Levitt α-helix propensity

**Pipeline Steps:**

1. **Feature Scaling:** StandardScaler on training data

2. **Feature Selection via RandomForest:**
   - Train RF with 150 trees
   - Rank features by importance
   - Visualize top 15 features
   - **Top features:** hydrophobicity and TM propensity (biologically sensible!)

3. **Hyperparameter Optimization:**
   ```
   GridSearchCV parameters:
   - Kernels: linear, RBF (γ: scale, auto, 0.01, 0.1), polynomial (degree: 2, 3)
   - C values: [0.1, 1, 10, 100]
   - 5-fold cross-validation
   - Scoring metric: F1 (for imbalanced data)
   - Total: 36 parameter combinations tested
   ```

4. **Best Parameters Found:**
   - Kernel: RBF
   - C: 10
   - γ: 0.1
   - Probability calibration enabled

**Test Set Performance:**

| Class | Precision | Recall | F1-Score |
|-------|-----------|--------|----------|
| Negative | 0.98 | 0.98 | 0.98 |
| Positive | 0.87 | 0.85 | 0.86 |
| **Overall** | | | **0.97 (accuracy)** |
| **MCC** | | | **0.847** |

**Transmembrane Protein Analysis:**
- Total TM negative proteins in test set: 165
- False Positives from TM subset: 21 / 165
- **FPR on TM proteins: 0.127** (vs. 0.206 for Von Heijne)
- Shows SVM better distinguishes SP from TM N-termini

**Outputs:**
- `svm_augmented_feature_importances.png`
- `svm_confusion_matrix.png`
- `svm_roc_curve.png` (AUC = 0.98)
- `svm_precision_recall_curve.png`
- `svm_signal_peptide_model.joblib` - Trained model
- `feature_scaler.joblib` - Feature scaler

---

### 6. **deep_learning/** - Phase F: Deep Learning Approach

**Purpose:** Leverage neural networks to capture sequence patterns

**Scripts:**
- `deep_learning.py` - CNN + BiLSTM architecture
- `vonHeijne.py`, `SVM.py` - Comparison models
- `complete_comparison.py` - Final benchmarking

**Model Architecture:**

```
Input: Protein sequence (variable length)
  ↓
One-Hot Encoding: 20 × 70 tensor (padded/truncated to 70 residues)
  ↓
CNN Layer:
  - 64 filters, kernel size 17
  - Captures local amino acid motifs
  - ReLU activation + Batch Normalization
  ↓
Bi-LSTM Layer:
  - 2 layers, 128 hidden units each
  - Learns long-range dependencies in both directions
  - Dropout (0.5) for regularization
  ↓
Fully Connected Layers:
  - Flatten LSTM output (take last time step)
  - 256 → 128 → 1 (binary output)
  - Dropout between layers
  ↓
Output: Sigmoid activation (probability)
```

**Training Details:**

- **Data Split:** 80% train, 10% validation, 10% test
- **Loss Function:** Binary Crossentropy with class weights (handles imbalance)
- **Optimizer:** Adam (default LR)
- **Early Stopping:** Patience=15 epochs, monitored on validation MCC
- **Batch Size:** 32
- **Epochs:** Up to 100 (typically stops earlier)

**Test Set Performance:**

| Class | Precision | Recall | F1-Score |
|-------|-----------|--------|----------|
| Negative | 0.99 | 0.99 | 0.99 |
| Positive | 0.94 | 0.93 | 0.93 |
| **Overall** | | | **0.99 (accuracy)** |
| **MCC** | | | **0.92** |

**Transmembrane Protein Analysis:**
- False Positives from TM subset: **6 / 165**
- **FPR on TM proteins: 0.036** ← **BEST PERFORMANCE**
- Significantly outperforms both SVM and Von Heijne!

**Outputs:**
- `best_signal_peptide_model.pth` - Best model checkpoint
- Confusion matrix, ROC curve, PR curve plots
- Model achieves highest accuracy on both standard and TM protein subsets


---

## Running the Pipeline

### Prerequisites

```bash
# Required packages
pip install requests pandas numpy scikit-learn torch matplotlib seaborn joblib tabulate
```

### Step-by-Step Execution

#### 1. Data Collection (Phase A)
```bash
cd data-collection
python dataset_gathering.py
# Outputs: positive_dataset_sp_cleavage.{tsv,fasta}, negative_dataset.{tsv,fasta}
```

#### 2. Data Preprocessing (Phase B)
```bash
# Note: MMseqs2 clustering requires VM setup (see mmseq_results/README.md)
# Assuming cluster results are available:
cd mmseq_results
python PartB_mmseq_preprocess.py
# Outputs: merged_dataset_with_seqs.tsv
```

#### 3. Data Analysis (Phase C)
```bash
cd data-analysis
python PartC_data_analysis.py
# Outputs: plots/
```

#### 4. Von Heijne Model (Phase D)
```bash
cd model-creation
python model_creation.py
# Outputs: cross_validation_results.txt, plots/
```

#### 5. SVM Model (Phase E)
```bash
cd SVM
python SVM.py
# Outputs: svm_signal_peptide_model.joblib, feature_scaler.joblib, plots/
```

#### 6. Deep Learning Model (Phase F)
```bash
cd deep_learning
python deep_learning.py
# Outputs: best_signal_peptide_model.pth, plots/
```

#### 7. Model Comparison
```bash
cd deep_learning
python complete_comparison.py
# Side-by-side evaluation of all three models
```

---

## Dependencies

### Core Libraries
- **Python 3.7+**
- **pandas**: Data manipulation and analysis
- **NumPy**: Numerical computing
- **scikit-learn**: Machine learning (SVM, RandomForest, metrics)
- **PyTorch**: Deep learning framework
- **Matplotlib & Seaborn**: Visualization
- **Requests**: HTTP API client

### Optional
- **MMseqs2**: Sequence clustering (requires separate installation)
- **WebLogo**: Sequence logo visualization

---
