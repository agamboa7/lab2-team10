# Data Analysis 

This part of the project was carried out in order to analyze and compare the distribution of the training and benchmarking data sets, which were obtained in the previous phase of the project and can be found in the `mmseq_results` directory of this repository. The datasets were compared in respect to the: protein length, signal peptide length, amino acid composition of signal peptides, taxonomic classification (kingdom, species) and signal peptide cleavage sites.

The programming tools that were used include Pandas (to create the data frames for further analysis) and Seaborn (to produce graphs that facilitate statistical data visualization). 

The programming tools that were used include Pandas (to create the data frames for further analysis) and Seaborn (to produce graphs that facilitate statistical data visualization). 


## Distribution of Protein Lenghts
To represent the distribution of protein lengths in the positive vs the negative sets in the testing and training datasets we used a density plot in the Seaborn library, plotting protein length on x-axis against the density on the y-axis. A kernel density estimate (KDE) plot is a method for visualizing the distribution of observations in a dataset, analogous to a histogram. KDE represents the data using a continuous probability density curve in one or more dimensions.
...

## Distribution of SP Lengths
To represent the distribution of protein lengths in signal peptides we used a histogram in the Seaborn library, plotting the cleavage site position on the x-axis against the density on the y-axis. A histogram is a classic visualization tool that represents the distribution of one or more variables by counting the number of observations that fall within discrete bins.
...

## Amino-acid Composition
The amino-acid composition of the signal peptides was compared against the amino acid composition of Swiss-prot, used as the background distribution. The amino acid composition of SwissProt is available at https://web.expasy.org/docs/relnotes/relstat.html
...

## Taxonomic Classification

We produced two pie charts (one for training set and one for test set) in order to compare the composition percentages at kingdom level in both sets. The resulting percentages are shown in the following table, showing that the distributions of kingdoms in both sets are consistent:

| Kingdom  | Train Set (%) | Test Set (%) |
|----------|---------------|--------------|
| Metazoa  | 55.4          | 55.9         |
| Fungi    | 25.7          | 25.5         |
| Plants   | 17.0          | 16.6         |
| Other    | 2.0           | 1.9          |


## Sequence Logos

...

