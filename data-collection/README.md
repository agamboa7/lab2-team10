# Data Collection

From the UniprotKB we retrieved the preliminary positive and negative datasets. In order to retrieve the url for the positive set we used the following filters in the advanced search: 

* Select only eukaryotic proteins
* Filter-out sequences shorter than 40 residues
* Filter-out unreviewed proteins
* Select on protein with experimental SP evidence
* No fragments

Query: (fragment:false) AND (taxonomy_id:2759) AND (length:[40 TO *]) AND (reviewed:true) AND (existence:1) AND (ft_signal_exp:*)
Two extra filters (proteins with SP shorter than 14 residues and cleaved) were implemented in our script: "script_name", function extract_fields.

In order to retrieve the url for the negative set we used the following filters in the advanced search: 

* No fragments
* Select only eukaryotic proteins
* Filter-out sequences shorter than 40 residues
* Filter-out sequences having SP (any evidence)
* Select only proteins experimentally verified to be localized into: cytosol, nucleus, mitochondrion, plastid, peroxisome, cell membrane

Query: (fragment:false) AND (reviewed:true) AND (existence:1) AND (taxonomy_id:2759) AND (length:[40 TO *]) AND ((cc_scl_term_exp:SL-0091) OR (cc_scl_term_exp:SL-0191) OR (cc_scl_term_exp:SL-0173) OR (cc_scl_term_exp:SL-0209) OR (cc_scl_term_exp:SL-0204) OR (cc_scl_term_exp:SL-0039)) NOT (ft_signal:*)

The script "script_name" takes as input the url for negative/positive dataset and gives as output the data stored in the tsv file (tsv file contains protein UniProt accession
number, the organism name, the Eukaryotic kingdom, the protein length, the position of the signal peptide cleavage site) and sequences stored in a FASTA file.
