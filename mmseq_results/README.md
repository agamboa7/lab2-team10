# README

To get these results we performed the following steps:

### **1. Perform MMseqs2 Clustering (on VM)**

- First, we performed MMseqs2 clustering inside the virtual machines (VMs) provided by the professor.
- As input, we used the `positive_dataset_sp_cleavage.fasta` and `negative_dataset.fasta` files.
- The specific MMseqs2 command used was:
  ```bash
  # [Insert your MMseqs2 clustering script/command here, e.g.,]
  # mmseqs easy-cluster positive_dataset_sp_cleavage.fasta pos_cluster tmp --min-seq-id 0.9 --cluster-mode 0
  # mmseqs easy-cluster negative_dataset.fasta neg_cluster tmp --min-seq-id 0.9 --cluster-mode 0
  ```
- After MMseqs2 finished, we retrieved the resulting `.tsv` files (`pos_cluster-results_cluster.tsv` and `neg_cluster-results_cluster.tsv`) from the VM to complete the rest of the analysis on Colab.

### **2. Select Representative Sequences**

- Using the clustered data (`pos_cluster-results_cluster.tsv` and `neg_cluster-results_cluster.tsv`), we identified the main "representative" sequence from each cluster.
- We then filtered the original datasets (`positive_dataset.tsv` and `negative_dataset.tsv`) to keep only these representative sequences.

### **3. Create Filtered Datasets**

- The newly filtered positive and negative datasets are saved as `filtered_pos_dataset.tsv` and `filtered_neg_dataset.tsv`.

### **4. Split Data into Training and Test Sets**

- The filtered data is split into a training set (80%) and a testing set (20%).
- A `label` is added: `1` for positive examples and `0` for negative examples.
- The final training and test sets (`train_df` and `test_df`) are created by combining the positive and negative splits.

### **5. Prepare for Cross-Validation**

- The training data (`train_df`) is divided into 5 "folds" for stratified cross-validation. This ensures that each fold has a similar ratio of positive to negative examples.

### **6. Merge and Add Sequences**

- The training and test sets are combined into a single file: `merged_dataset.tsv`.
- The actual protein sequences are read from the original FASTA files and added to this merged file.
- The final result, `merged_dataset_with_seqs.tsv`, contains all the prepared data, including sequences.

### **7. Generate a Summary**

- At the end, the script prints a summary table showing the number of proteins at each stage of the process (before and after clustering, and in the final train/test splits).


| Dataset                      |   Total Samples |   Negative |   Positive |
|:-----------------------------|----------------:|-----------:|-----------:|
| Merged (Train)               |            8021 |       7147 |        874 |
| Merged (Test)                |            2006 |       1787 |        219 |
| Negative (Before Clustering) |           20615 |      20615 |          0 |
| Positive (Before Clustering) |            2932 |          0 |       2932 |
| Negative (After Clustering)  |            8934 |       8934 |          0 |
| Positive (After Clustering)  |            1093 |          0 |       1093 |
