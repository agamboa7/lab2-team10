# Data Collection

From the UniprotKB we retrieved the preliminary positive and negative datasets. In order to retrieve the url for the positive set we used the following filters in the advanced search: 

* Select only eukaryotic proteins
* Filter-out sequences shorter than 40 residues
* Filter-out unreviewed proteins
* Select on protein with experimental SP evidence
* No fragments

Query: (fragment:false) AND (taxonomy_id:2759) AND (length:[40 TO *]) AND (reviewed:true) AND (existence:1) AND (ft_signal_exp:*)

Two extra filters (proteins with SP shorter than 14 residues and cleaved) were implemented in our script: "DataCollection_team10", function extract_fields.

In order to retrieve the url for the negative set we used the following filters in the advanced search: 

* No fragments
* Select only eukaryotic proteins
* Filter-out sequences shorter than 40 residues
* Filter-out sequences having SP (any evidence)
* Select only proteins experimentally verified to be localized into: cytosol, nucleus, mitochondrion, plastid, peroxisome, cell membrane

Query: (fragment:false) AND (reviewed:true) AND (existence:1) AND (taxonomy_id:2759) AND (length:[40 TO *]) AND ((cc_scl_term_exp:SL-0091) OR (cc_scl_term_exp:SL-0191) OR (cc_scl_term_exp:SL-0173) OR (cc_scl_term_exp:SL-0209) OR (cc_scl_term_exp:SL-0204) OR (cc_scl_term_exp:SL-0039)) NOT (ft_signal:*)

The script "DataCollection_team10" takes as input the url for negative/positive dataset and gives as output the data stored in the tsv file (tsv file contains protein UniProt accession number, the organism name, the Eukaryotic kingdom, the protein length, the position of the signal peptide cleavage site) and sequences stored in a FASTA file.

## Extra filtering on the script

Since there were two criteria that cannot be obtained through the advanced search tool in UniProtKB, the resulting entries (2949 proteins) were additionally filtered using the features on the json and retrieving the ones that were cleaved (identified by a characteristic empty space on the description of the signal peptide) and filtering out the ones with signal peptide shorter than 14 residues. The filtering was performed with the following function:

~~~
def extract_fields(entry):

  accession = entry["primaryAccession"]
  organism_name = entry["organism"]["scientificName"]
  lineage = entry["organism"].get("lineage", [])
  kingdom = get_eukaryotic_kingdom(lineage)
  protein_length = entry["sequence"]["length"]

  cleavage_site = 0 # Default value if not found
  # Find the signal peptide feature to get the cleavage site
  for f in entry.get("features", []):
    if f.get("type") == "Signal":
      if f["description"] != "":
        return None # Signal to filter this entry out
      if f["location"]["end"]["value"] == "?":
        return None # Signal to filter this entry out
      if not f:
        return None

      # The cleavage site is the end position of the signal peptide
      start = f["location"]["start"]["value"]
      end = f["location"]["end"]["value"]
      cleavage_len = end - start + 1
      if cleavage_len < 14:
        return None # Signal to filter this entry out
      cleavage_site = cleavage_len
      break # Stop after finding the first one

  return (
      accession,
      organism_name,
      kingdom,
      protein_length,
      cleavage_site
  )
  ~~~

...
