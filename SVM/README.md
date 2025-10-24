# SVM with Feature Selection

This script builds a pipeline to classify proteins based on the presence of a signal peptide (SP), focusing on the implementation and optimization of a Support Vector Machine (SVM).


## Table of Contents
- [Script Overview](#script-overview)
- [Pipeline Steps](#pipeline-steps)
- [Generated Files](#generated-files)
- [Results](#key-results-from-this-session)

## Overview

The goal of the `SVM.py` script is to apply advanced machine learning techniques to the signal peptide classification problem. The script takes the pre-processed dataset (`merged_dataset_with_seqs.tsv`) and executes a pipeline that includes feature engineering, feature selection, hyperparameter tuning, model evaluation, and final model persistence. This represents the culmination of the modeling phase of the course project. 

Afterwards, the `model_comparison.py` script compares the performance of the newly created SVM model with the previous von Heijne model. Spesifically, we compare their performance on classifying transmembrane proteins.

## Pipeline Steps

The script automates the following sequence of tasks:

1.  **Data Loading:** Loads the training and testing data splits from the course project's dataset file.
2.  **Feature Engineering:** Creates a rich set of biochemical features from the raw N-terminal protein sequences, based on established physicochemical properties of amino acids.
3.  **Feature Importance Ranking:** A RandomForestClassifier is trained to evaluate the predictive power of each engineered feature.
4.  **Feature Filtering:** The dataset is reduced to include only the **top 15 most important features** as identified by the RandomForest.
5.  **SVM Hyperparameter Optimization:** `GridSearchCV` is used to systematically find the best hyperparameters for an SVM classifier (testing `linear`, `rbf`, and `poly` kernels) using the pre-defined cross-validation folds.
6.  **Final Model Evaluation:** The best-performing SVM model is evaluated on the hold-out test set to assess its generalization performance.
7.  **Model and Scaler Export:** The final, trained SVM model and the feature scaler are saved to disk using `joblib`, making them available for future use without retraining.


## Generated Files

Executing this script will produce the following output files in the working directory:

*   `svm_augmented_feature_importances.png`: A bar plot showing the importance of each feature, as determined by the RandomForest.
*   `svm_confusion_matrix.png`: A confusion matrix visualizing the performance of the final model on the test data.
*   `svm_roc_curve.png`: The Receiver Operating Characteristic (ROC) curve for the final model, including the Area Under the Curve (AUC) score.
*   `svm_signal_peptide_model.joblib`: The final, trained SVM classifier object.
*   `feature_scaler.joblib`: The `StandardScaler` object fitted to the training data.

## Key Results from this Session

| Classification | precision | recall | f1-score | support |
| :--- | :--- | :--- | :--- | :--- |
| **Negative (0)** | 0.98 | 0.98 | 0.98 | 1787 |
| **Positive (1)** | 0.87 | 0.85 | 0.86 | 219 |
| | | | | |
| **accuracy** | | | 0.97 | 2006 |
| **macro avg** | 0.93 | 0.92 | 0.92 | 2006 |
| **weighted avg** | 0.97 | 0.97 | 0.97 | 2006 |

The following plots summarize the main outcomes of this implementation step: the ranking of features, and the performance of the final, optimized SVM model on the test set.

<img src="Feature Importences from Random Forest.png" alt="Feature Importences from Random Forest" width="60%">

<img src="Confusion Matrix on Test Set.png" alt="COnfusion Matrix" width="60%">

<img src="Roc Curve.png" alt="Roc Curve" width="60%">


### Model Performance Comparison

#### 1. False Positive Rate (FPR) on Transmembrane Proteins

A key challenge is distinguishing signal peptides from the N-terminal anchor helices of transmembrane (TM) proteins, as both are hydrophobic. A lower FPR on this subset indicates a better model.

| Metric | SVM Model | Von Heijne Model |
| :--- | :--- | :--- |
| **Standard FPR (Overall)** | **`0.0151`** | `0.0330` |
| **FPR on TM Proteins** | **`0.1273`** | `0.2061` |
| False Positives from TM Subset | **`21`** / 165 | `34` / 165 |

> **Conclusion**: The **SVM model is significantly better** at correctly identifying transmembrane N-termini as negative. Its specialized FPR is much lower, suggesting it learned features that can more effectively distinguish between TM anchors and true signal peptides.

#### 2. Analysis of Misclassified Sequences

This analysis examines the average physicochemical properties of sequences that were misclassified by each model. The "h-region hydrophobicity" is the most critical feature for signal peptide recognition.

##### SVM Model: Error Analysis

| Prediction Outcome | Avg. n-region Charge | Avg. h-region Hydrophobicity | Count |
| :--- | :---: | :---: | :---: |
| **True Positives** | 0.545 | `2.183` | 187 |
| **False Negatives** | 0.219 | `0.330` | 32 |
| **True Negatives** | -0.005 | `-0.478`| 1760 |
| **False Positives** | 0.407 | `2.044` | 27 |

*   **Key Insight (False Positives)**: The SVM's false positives have an extremely high hydrophobicity (`2.044`), which is nearly identical to that of true positives (`2.183`). This shows the model is primarily confused by sequences that strongly mimic a signal peptide's hydrophobic core.
*   **Key Insight (False Negatives)**: The model misses true signal peptides that have unusually low hydrophobicity (`0.330`), making them atypical.

##### Von Heijne Model: Error Analysis

| Prediction Outcome | Avg. n-region Charge | Avg. h-region Hydrophobicity | Count |
| :--- | :---: | :---: | :---: |
| **True Positives** | 0.481 | `1.969` | 154 |
| **False Negatives** | 0.538 | `1.778` | 65 |
| **True Negatives** | 0.000 | `-0.459`| 1728 |
| **False Positives** | 0.034 | `0.120` | 59 |

*   **Key Insight (False Positives)**: The Von Heijne model's false positives have only a slightly elevated hydrophobicity (`0.120`) compared to true negatives. This suggests its errors are caused by more subtle sequence patterns rather than just high hydrophobicity.
*   **Key Insight (False Negatives)**: The model struggles with a larger number of signal peptides (`65`) that have a hydrophobicity (`1.778`) only slightly lower than the average true positive.


