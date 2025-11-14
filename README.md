Of course. Here is the updated README with all deep learning-related content removed:

# lab2-team10
This repository contains material relevant to the Laboratory of Bioinformatics 2 course, which is part of the Bioinformatics Master's Degree at the University of Bologna.

**Team Members:**
- Andrea Arriola Gamboa
- Kristina Djikic
- Başak Akkoyun
- Betül Yalçın
- Deniz Ertuğrul

# Signal Peptide Detection: Machine Learning Pipeline

This repository contains a comprehensive machine learning pipeline for predicting signal peptides (SP) in eukaryotic protein sequences. The project compares multiple approaches including the classic von Heijne algorithm and Support Vector Machines (SVM).

---

## Table of Contents

1. [Project Overview](#project-overview)
2. [Repository Structure](#repository-structure)
3. [Pipeline Architecture](#pipeline-architecture)
4. [Key Results](#key-results)
5. [Detailed Module Descriptions](#detailed-module-descriptions)
6. [Dependencies](#dependencies)

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

## Pipeline Architecture

```mermaid
graph LR
    A[" Data Collection"] --> B[" Preprocessing"]
    B --> C["Analysis"]
    C --> D["von Heijne Method"]
    C --> E["Support Vector Machine"]
    D --> G["Model Evaluation"]
    E --> G
    G --> H["Deployment"]
````

---

## Repository Structure

```
📦 lab2-team10/
│
├── 📁 data-collection/        # Phase A — UniProt API integration
├── 📁 mmseq_results/          # Phase B — Clustering & redundancy reduction
├── 📁 data-analysis/          # Phase C — Exploratory analysis & visualization
├── 📁 model-creation/         # Phase D — Von Heijne statistical model
└── 📁 SVM/                    # Phase E — Support Vector Machine
```

---

## Key Results

### Model Comparison Summary

| Metric | Von Heijne | SVM |
|--------|-----------|-----|
| **Overall Accuracy** | 0.94 | **0.97** |
| **F1 Score (Positive)** | 0.71 | **0.86** |
| **Precision (Positive)** | 0.72 | **0.87** |
| **Recall (Positive)** | 0.70 | **0.85** |
| **MCC** | 0.694 | **0.847** |
| **FPR on TM Proteins** | 0.206 | **0.127** |

### Key Findings

1. **SVM is the top performer:** It achieves the highest accuracy (97%) and best handles transmembrane proteins (FPR: 12.7%).

2. **Von Heijne is a competitive baseline:** The traditional method achieves 94% accuracy and serves as a useful reference.

3. **TM protein challenge:** Both models struggle more with TM proteins (higher false positive rate) due to similar N-terminal hydrophobicity. SVM mitigates this most effectively.

4. **Feature importance:** Hydrophobicity and TM propensity are the top predictive features for the SVM model, confirming the biological understanding of SP structure.

---

## Detailed Module Descriptions

### Phase A Data Collection

<details>
<summary><b>data-collection/</b> — Fetch and preprocess raw protein data from UniProtKB</summary>

**Script:** `dataset_gathering.py`

#### Key Features:

* Robust API client with retry logic for server errors
* Pagination handling for large search results (~3000 positive, ~20,600 negative)
* Advanced filtering for quality control
* Dual format output: TSV (metadata) + FASTA (sequences)

#### UniProt Search Queries:

**Positive:**

```
(fragment:false) AND (taxonomy_id:2759) AND (length:[40 TO *]) 
AND (reviewed:true) AND (existence:1) AND (ft_signal_exp:*)
```

**Negative:**

```
Secreted proteins with cytosolic/nuclear/organellar localization,
excluding those with signal peptides
```

#### Kingdom Distribution:

| Kingdom      | Metazoa | Plants | Fungi | Other |
| ------------ | ------- | ------ | ----- | ----- |
| **Positive** | 82.5%   | 10.6%  | 5.6%  | 1.2%  |
| **Negative** | 60.2%   | 19.9%  | 18.1% | 1.7%  |

#### Outputs:

* `positive_dataset_sp_cleavage.tsv` / `.fasta`
* `negative_dataset.tsv` / `.fasta`

</details>

---

### Phase B Data Preprocessing & Clustering

<details>
<summary><b>mmseq_results/</b> — Remove redundancy and create train/test splits</summary>

**Purpose:** Eliminate sequence redundancy and prepare data for model training

#### Key Steps:

**1. MMseqs2 Clustering** (on VM):

```bash
mmseqs easy-cluster input.fa cluster-results tmp \
  --min-seq-id 0.3 -c 0.4 --cov-mode 0
```

* 30% minimum sequence identity
* Reduces 23,547 → 10,027 sequences

**2. Representative Selection:**

* Extract cluster representatives (highest similarity to all members)
* Final: **1,093 positive** + **8,934 negative**

**3. Train/Test Split:**

* 80% training / 20% test (stratified)
* Ensures class balance in both sets

**4. 5-Fold Cross-Validation:**

* Training set divided into 5 balanced folds
* Used for hyperparameter tuning

**5. Sequence Integration:**

* Merge FASTA with metadata
* Output: `merged_dataset_with_seqs.tsv`

#### Dataset Statistics:

| Set       | Negative  | Positive  | Total      |
| --------- | --------- | --------- | ---------- |
| Train     | 7,147     | 874       | 8,021      |
| Test      | 1,787     | 219       | 2,006      |
| **Total** | **8,934** | **1,093** | **10,027** |

</details>

---

### Phase C Exploratory Data Analysis

<details>
<summary><b>data-analysis/</b> — Characterize dataset properties</summary>

**Script:** `PartC_data_analysis.py`

#### Analyses Performed:

1. **Protein Length Distributions**

   * Positive vs. negative sequences
   * Median positive: ~300 aa | Median negative: ~450 aa
   * Negative sequences are significantly longer

2. **Signal Peptide Length Analysis**

   * Cleavage site positions
   * Most SPs: 14–30 residues (typical range)

3. **Amino Acid Composition**

   * Compare SP sequences vs. SwissProt background
   * **Enriched:** Small hydrophobic (Ala, Val, Leu)
   * **Depleted:** Charged residues (Arg, Asp, Glu)

4. **Taxonomic Distribution**

   * Kingdom-level breakdown
   * Consistent train/test distributions

5. **Cleavage Site Motifs**

   * Extract [-13, +2] windows
   * Generate sequence logos
   * Reveal conserved patterns

#### Key Finding:

Train and test sets have **consistent distributions** across all metrics.

</details>

---

### Phase D Von Heijne Algorithm

<details>
<summary><b>model-creation/</b> — Traditional statistical approach</summary>

**Script:** `model_creation.py`

#### Algorithm Overview:

**1. Extract Training Motifs**

* From positive training proteins: [-13, +2] windows at cleavage sites

**2. Build PSWM (Position-Specific Probability Matrix)**
$$PSPM[aa][pos] = \frac{count[aa] + pseudocount}{total + background}$$

**3. Compute Weight Matrix**
$$PSWM[aa][pos] = \log\left(\frac{PSPM[aa][pos]}{background_freq[aa]}\right)$$

**4. Prediction**

* Scan N-terminus (first 90 residues)
* Score all windows of length 16
* Maximum score = protein's SP score

**5. Threshold Optimization**

* F1-score maximization via Precision-Recall curve
* Optimal threshold ≈ 6.73

#### Cross-Validation Results:

| Metric        | Value         |
| ------------- | ------------- |
| **F1 Score**  | 0.724 ± 0.016 |
| **Precision** | 0.750 ± 0.014 |
| **Recall**    | 0.706 ± 0.035 |
| **MCC**       | 0.694 ± 0.016 |

#### Outputs:

* `averaged_precision_recall_curve.png`
* `pswm_heatmap.png`

</details>

---

### Phase E Support Vector Machine

<details>
<summary><b>SVM/</b> — Advanced ML with feature engineering</summary>

**Scripts:** `SVM.py` | `vonHeijne.py` (comparison) | `model_comparison.py`

#### Feature Engineering (26 Total):

**Amino Acid Composition (20 features):**

* Frequency of each standard amino acid in N-terminal region (22 residues)

**Physicochemical Properties (6 features):**

* n-region charge (first 5 residues)
* h-region hydrophobicity (residues 5–16)
* h-region TM propensity (residues 5–16)
* Total N-terminal hydrophobicity
* Total N-terminal TM propensity
* Total N-terminal α-helix propensity

#### Pipeline:

**1. Feature Scaling:**

```python
StandardScaler(fit_params=training_data)
```

**2. Feature Selection via Random Forest:**

* Train RF (150 trees)
* Rank by importance
* Top features: **Hydrophobicity & TM propensity**

**3. Hyperparameter Optimization:**

```python
GridSearchCV(
    kernels=['linear', 'RBF', 'polynomial'],
    C_values=[0.1, 1, 10, 100],
    gamma=['scale', 'auto', 0.01, 0.1],
    cv=5,
    scoring='f1'
)
```

**4. Best Parameters:**

* Kernel: **RBF**
* C: **10**
* γ: **0.1**

#### Test Performance:

| Class             | Precision | Recall | F1-Score  |
| ----------------- | --------- | ------ | --------- |
| Negative          | 0.98      | 0.98   | 0.98      |
| Positive          | 0.87      | 0.85   | 0.86      |
| **Overall (Acc)** |           |        | **0.97**  |
| **MCC**           |           |        | **0.847** |

#### TM Protein Analysis:

* TM negatives in test: 165
* False positives: 21 / 165
* **FPR: 0.127** ← Better than Von Heijne (0.206)

#### Outputs:

* `svm_augmented_feature_importances.png`
* `svm_confusion_matrix.png`
* `svm_roc_curve.png` (AUC = 0.98)
* `svm_precision_recall_curve.png`
* `svm_signal_peptide_model.joblib`
* `feature_scaler.joblib`

</details>

---

## Dependencies

### Core Libraries
- **Python 3.7+**
- **pandas**: Data manipulation and analysis
- **NumPy**: Numerical computing
- **scikit-learn**: Machine learning (SVM, RandomForest, metrics)
- **Matplotlib & Seaborn**: Visualization
- **Requests**: HTTP API client

### Optional
- **MMseqs2**: Sequence clustering (requires separate installation)
- **WebLogo**: Sequence logo visualization

---

### 🌟 Made with ❤️ by Team 10

**University of Bologna** — Bioinformatics Master's Program
