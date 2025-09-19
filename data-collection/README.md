# Data Collection

From the UniprotKB the preliminary positive and negative datasets were retrieved. In order to obtain the URL for the positive set, the following criteria were applied in the advanced search: 

* Select only eukaryotic proteins
* Filter-out sequences shorter than 40 residues
* Filter-out unreviewed proteins
* Select on protein with experimental SP evidence
* No fragments

Query: (fragment:false) AND (taxonomy_id:2759) AND (length:[40 TO \*]) AND (reviewed:true) AND (existence:1) AND (ft_signal_exp:\*)

Two extra filters (proteins with SP shorter than 14 residues and cleaved) were implemented in the script: `dataset_gathering.py`, function `extract_fields_positive`.

In order to retrieve the URL for the negative set, the following criteria was considered in the advanced search: 

* No fragments
* Select only eukaryotic proteins
* Filter-out sequences shorter than 40 residues
* Filter-out sequences having SP (any evidence)
* Select only proteins experimentally verified to be localized into: cytosol, nucleus, mitochondrion, plastid, peroxisome, cell membrane

Query: (fragment:false) AND (reviewed:true) AND (existence:1) AND (taxonomy_id:2759) AND (length:[40 TO \*]) AND ((cc_scl_term_exp:SL-0091) OR (cc_scl_term_exp:SL-0191) OR (cc_scl_term_exp:SL-0173) OR (cc_scl_term_exp:SL-0209) OR (cc_scl_term_exp:SL-0204) OR (cc_scl_term_exp:SL-0039)) NOT (ft_signal:\*)

The script `dataset_gathering.py` takes as input the URL for negative/positive datasets and gives as output the data stored in the TSV file (which contains protein UniProt accession number, organism name, Eukaryotic kingdom, protein length, and position of the signal peptide cleavage site) and sequences stored in a FASTA file.

## Extra filtering on the script

Since there were two criteria that could not be obtained through the advanced search tool in UniProtKB, the resulting entries (2,949 proteins) were additionally filtered using the features on the JSON and retrieving the ones that were cleaved (identified by a characteristic empty space on the description of the signal peptide) and filtering out the ones with signal peptide shorter than 14 residues. The filtering was performed with the following function:

~~~
def extract_fields_positive(entry):
    """
    Extracts the required fields for the positive dataset.
    
    Filters out entries that don't meet the criteria for a good positive example.
    """
    accession = entry["primaryAccession"]
    organism_name = entry["organism"]["scientificName"]
    lineage = entry["organism"].get("lineage", [])
    kingdom = get_eukaryotic_kingdom(lineage)
    protein_length = entry["sequence"]["length"]
    sequence = entry["sequence"]["value"]

    cleavage_site = 0  # Default value if not found
    # Find the signal peptide feature to get the cleavage site
    for f in entry.get("features", []):
        if f.get("type") == "Signal":
            # Filter out entries with descriptive text, missing end location, or short signal peptides
            if f.get("description", "") != "" or f["location"]["end"]["value"] == "?":
                return None
            
            start = f["location"]["start"]["value"]
            end = f["location"]["end"]["value"]
            cleavage_len = end - start + 1
            if cleavage_len < 14:
                return None  # Filter out entries with a short signal peptide
            cleavage_site = cleavage_len
            break  # Stop after finding the first one

    return (
        accession,
        organism_name,
        kingdom,
        protein_length,
        cleavage_site,
        sequence
    )
  ~~~

> Complete script available on this directory under the name `dataset_gathering.py`

The python script successfully generated the required datasets for the project, which were also uploaded on the `data-collection` directory of this repository. The generated files are the following:

* positive_dataset_sp_cleavage.tsv
* positive_dataset_sp_cleavage.fasta
* negative_dataset.tsv
* negative_dataset.fasta

## To run the script on your own:

Make sure to have python3, download the `dataset_gathering.py` script, and run the following on the command line:

```bash
python3 dataset_gathering.py
``` 
- `output`: 4 output files (2 TSV files and 2 FASTA files, each for positive and negative datasets)


## Final Dataset Summary Table

| Dataset         | Number of Proteins | Metazoa | Plants | Fungi | Other | Proteins with N-terminal TMH |
|-----------------|--------------------|---------|--------|-------|-------|------------------------------|
| Positive        | 2932                 | 2420      | 311     | 165    | 36    | N/A                          |
| Negative        | 20615                 | 12419      | 4111     | 3727    | 358    | 2477                           |

