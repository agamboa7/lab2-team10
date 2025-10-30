# Final Model Comparison: SVM vs. Deep Learning vs. Von Heijne

## 1. Data Loading
- Data loaded from `merged_dataset_with_seqs.tsv`.
- **Training set size:** 8021
- **Test set size:** 2006

## 2. Model Predictions

## 3. Saving Combined Results
- Extended test set with all model predictions saved to `test_set_with_all_predictions.tsv`.

## 4. Comparative Performance Analysis

### SVM Model: Classification Report
| Metric          | Precision | Recall | F1-Score | Support |
|-----------------|-----------|--------|----------|---------|
| **Negative**       | 0.98      | 0.98   | 0.98     | 1787.0      |
| **Positive**       | 0.87      | 0.85   | 0.86     | 219.0      |
|-----------------|-----------|--------|----------|---------|
| **Accuracy**      |           |        | 0.97     | 2006.0      |
| **Macro Avg**     | 0.93      | 0.92   | 0.92     | 2006.0      |
| **Weighted Avg**  | 0.97      | 0.97   | 0.97     | 2006.0      |

### DL Model: Classification Report
| Metric          | Precision | Recall | F1-Score | Support |
|-----------------|-----------|--------|----------|---------|
| **Negative**       | 0.99      | 0.99   | 0.99     | 1787.0      |
| **Positive**       | 0.94      | 0.93   | 0.93     | 219.0      |
|-----------------|-----------|--------|----------|---------|
| **Accuracy**      |           |        | 0.99     | 2006.0      |
| **Macro Avg**     | 0.96      | 0.96   | 0.96     | 2006.0      |
| **Weighted Avg**  | 0.99      | 0.99   | 0.99     | 2006.0      |

### Von_Heijne Model: Classification Report
| Metric          | Precision | Recall | F1-Score | Support |
|-----------------|-----------|--------|----------|---------|
| **Negative**       | 0.96      | 0.97   | 0.97     | 1787.0      |
| **Positive**       | 0.72      | 0.70   | 0.71     | 219.0      |
|-----------------|-----------|--------|----------|---------|
| **Accuracy**      |           |        | 0.94     | 2006.0      |
| **Macro Avg**     | 0.84      | 0.84   | 0.84     | 2006.0      |
| **Weighted Avg**  | 0.94      | 0.94   | 0.94     | 2006.0      |

## 5. Specialized Error Analysis

### SVM Model: Transmembrane Protein Analysis
- **Total Negative TM Proteins:** 165
- **False Positives from TM subset:** 21
- **FPR for Transmembrane proteins:** `0.1273`

### DL Model: Transmembrane Protein Analysis
- **Total Negative TM Proteins:** 165
- **False Positives from TM subset:** 6
- **FPR for Transmembrane proteins:** `0.0364`

### Von_Heijne Model: Transmembrane Protein Analysis
- **Total Negative TM Proteins:** 165
- **False Positives from TM subset:** 34
- **FPR for Transmembrane proteins:** `0.2061`

### SVM Model: Misclassification Analysis
- **Average H-region Hydrophobicity of True Positives:** `2.183`
- **Average H-region Hydrophobicity of FALSE NEGATIVES:** `0.330`
- **Average H-region Hydrophobicity of FALSE POSITIVES:** `2.044`

**Interpretation:**
- FN hydrophobicity is lower than TP, suggesting the model misses less canonical (less hydrophobic) signal peptides.
- FP hydrophobicity > 0 suggests the model is confused by hydrophobic N-termini (like TM anchors).

### DL Model: Misclassification Analysis
- **Average H-region Hydrophobicity of True Positives:** `1.999`
- **Average H-region Hydrophobicity of FALSE NEGATIVES:** `0.731`
- **Average H-region Hydrophobicity of FALSE POSITIVES:** `1.297`

**Interpretation:**
- FN hydrophobicity is lower than TP, suggesting the model misses less canonical (less hydrophobic) signal peptides.
- FP hydrophobicity > 0 suggests the model is confused by hydrophobic N-termini (like TM anchors).

### Von_Heijne Model: Misclassification Analysis
- **Average H-region Hydrophobicity of True Positives:** `1.969`
- **Average H-region Hydrophobicity of FALSE NEGATIVES:** `1.778`
- **Average H-region Hydrophobicity of FALSE POSITIVES:** `0.120`

**Interpretation:**
- FN hydrophobicity is lower than TP, suggesting the model misses less canonical (less hydrophobic) signal peptides.
- FP hydrophobicity > 0 suggests the model is confused by hydrophobic N-termini (like TM anchors).
