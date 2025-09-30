# Data Analysis 

This part of the project was carried out in order to analyze and compare the distribution of the training and benchmarking data sets, which were obtained in the previous phase of the project and can be found in the `mmseq_results` directory of this repository. The datasets were compared in respect to the: protein length, signal peptide length, amino acid composition of signal peptides, taxonomic classification (kingdom, species) and signal peptide cleavage sites.


The programming tools that were used include Pandas (to create the data frames for further analysis) and Seaborn (to produce graphs that facilitate statistical data visualization). The python script that was used throughout this step of the project is named `PartC_data_analysis.py` and is available in this repository.


## Distribution of Protein Lenghts
To represent the distribution of protein lengths in the positive vs the negative sets in the testing and training datasets we used a density plot in the Seaborn library, plotting protein length on x-axis against the density on the y-axis. A kernel density estimate (KDE) plot is a method for visualizing the distribution of observations in a dataset, analogous to a histogram. KDE represents the data using a continuous probability density curve in one or more dimensions.

[Click here to access the density plot of protein lengths in the whole dataset](https://github.com/agamboa7/lab2-team10/blob/main/data-analysis/plots/Density_Plot_of_Protein_Lengths_In_The_Whole_Dataset.png)

[Click here to access the density_plot_of_protein_lengths_in_train_and_test_sets](https://github.com/agamboa7/lab2-team10/blob/main/data-analysis/plots/Density_Plot_of_Protein_Lengths_in_Train_and_Test_sets.png)

[Click here to access the density_plot_of_protein_lengths_in_train_set](https://github.com/agamboa7/lab2-team10/blob/main/data-analysis/plots/Density_Plot_of_Protein_Lengths_in_Train_Set.png)

[Click here to access the density_plot_of_protein_lengths_in_test_set](https://github.com/agamboa7/lab2-team10/blob/main/data-analysis/plots/Density_Plot_of_Protein_Lengths_in_Test_Set.png)

[Click here to access the boxplot of protein lengths by train vs test](https://github.com/agamboa7/lab2-team10/blob/main/data-analysis/plots/Boxplot_of_Protein_Lengths_by_Train_vs_Test.png)


## Distribution of SP Lengths
To represent the distribution of protein lengths in signal peptides we used a histogram in the Seaborn library, plotting the cleavage site position on the x-axis against the density on the y-axis. A histogram is a classic visualization tool that represents the distribution of one or more variables by counting the number of observations that fall within discrete bins.

[Click here to access probability plot of SP lengths of train vs test](https://github.com/agamboa7/lab2-team10/blob/main/data-analysis/plots/Probability_Plot_of_SP_lengths_of_Train_vs_Test.png)


## Amino-acid Composition
The amino-acid composition of the signal peptides was compared against the amino acid composition of Swiss-prot, used as the background distribution. The amino acid composition of SwissProt is available at https://web.expasy.org/docs/relnotes/relstat.html

[Click here to access plot of AA composition](https://github.com/agamboa7/lab2-team10/blob/main/data-analysis/plots/AA_Composition_SwissProt_Train_Test_grouped.png)


## Taxonomic Classification

We produced two pie charts (one for training set and one for test set) in order to compare the composition percentages at kingdom level in both sets. The resulting percentages are shown in the following table, showing that the distributions of kingdoms in both sets are consistent:

| Kingdom  | Train Set (%) | Test Set (%) |
|----------|---------------|--------------|
| Metazoa  | 55.4          | 55.9         |
| Fungi    | 25.7          | 25.5         |
| Plants   | 17.0          | 16.6         |
| Other    | 2.0           | 1.9          |

[Click here to access plot of taxonomy classification](https://github.com/agamboa7/lab2-team10/blob/main/data-analysis/plots/Taxonomy_Classification_Train_vs_Test.png)

[Click here to access the pie chart of kingdom classification](https://github.com/agamboa7/lab2-team10/blob/main/data-analysis/plots/Kingdom_Classification_Train_vs_Test.png)

[Click here to access pie chart of organism classification](https://github.com/agamboa7/lab2-team10/blob/main/data-analysis/plots/Organism_Classification_Train_vs_Test.png)


## Sequence Logos

To analyze conserved motifs around the cleavage site, we extracted subsequences spanning positions [-13, +2] relative to the cleavage site. These motifs were saved in a text file and used as input to WebLogo, where the height of each amino acid represents its frequency at that position. This highlights conserved features typical of signal peptide cleavage regions.


