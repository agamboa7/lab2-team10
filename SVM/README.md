# SVM with Feature Selection

This script builds a pipeline to classify proteins based on the presence of a signal peptide (SP), focusing on the implementation and optimization of a Support Vector Machine (SVM).


## Table of Contents
- [Script Overview](#script-overview)
- [Pipeline Steps](#pipeline-steps)
- [Generated Files](#generated-files)
- [Results](#key-results-from-this-session)

## Script Overview

The goal of this script is to apply advanced machine learning techniques to the signal peptide classification problem. The script takes the pre-processed dataset (`merged_dataset_with_seqs.tsv`) and executes a pipeline that includes feature engineering, feature selection, hyperparameter tuning, model evaluation, and final model persistence. This represents the culmination of the modeling phase of the course project.

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

<img src="Roc Curve.png" alt="Roc Curve" width="60%">
