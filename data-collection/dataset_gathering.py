"""
This script fetches and processes protein data from the Uniprot API to create
positive and negative datasets for signal peptide analysis.

The script performs the following steps:
1. Defines a session with retry logic for robust API calls.
2. Implements a pagination mechanism to handle large search results.
3. Fetches a positive dataset of proteins with experimental signal peptides.
4. Extracts specific fields (including sequence) and filters the data based on defined criteria.
5. Fetches a negative dataset of proteins that are secreted but lack a
   signal peptide.
6. Extracts relevant fields (including sequence) from the negative dataset, including
   information about N-terminal transmembrane helices.
7. Saves both datasets to separate TSV and FASTA files.
8. Reads the generated TSV files into pandas DataFrames for a final summary.
"""

import requests
import re
import json
import pandas as pd
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

# We use the Session object to allow retries in case of temporary service unavailability.
# This configures it to try up to 5 times if it encounters specific server errors.
retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
session = requests.Session()
session.mount("https://", HTTPAdapter(max_retries=retries))

def get_next_link(headers):
    """
    Extracts the URL for the next page of results from the API response headers.
    """
    if "Link" in headers:
        # The regular expression is used to extract the next link for pagination
        re_next_link = re.compile(r'<(.+)>; rel="next"')
        match = re_next_link.match(headers["Link"])
        if match:
            return match.group(1)
    return None

def get_batch(batch_url):
    """
    Retrieves data one batch (page) at a time from a given URL.
    """
    while batch_url:
        # Run the API call
        response = session.get(batch_url)
        # Will raise an error if an unsuccessful status code is obtained
        response.raise_for_status()
        # Get the total number of entries in the search
        total = response.headers.get("x-total-results", "0")
        # Yield the response and the total number of entries
        yield response, total
        # Get the link to the API call for the next data batch
        batch_url = get_next_link(response.headers)

def get_eukaryotic_kingdom(lineage):
    """
    Categorizes an organism into Metazoa, Fungi, Plants, or Other based on its lineage.
    """
    if "Metazoa" in lineage:
        return "Metazoa"
    if "Fungi" in lineage:
        return "Fungi"
    if "Viridiplantae" in lineage:  # Viridiplantae is the taxonomic group for green plants
        return "Plants"
    return "Other"

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

def extract_fields_negative(entry):
    """
    Extracts the required fields for the negative dataset.
    """
    accession = entry["primaryAccession"]
    organism_name = entry["organism"]["scientificName"]
    lineage = entry["organism"].get("lineage", [])
    kingdom = get_eukaryotic_kingdom(lineage)
    protein_length = entry["sequence"]["length"]
    sequence = entry["sequence"]["value"]

    has_tmh_n_terminus = False  # Default value
    # Iterate through features to find a transmembrane helix (TMH)
    for f in entry.get("features", []):
        if f.get("type") == "Transmembrane":
            # Check if the TMH starts within the first 90 residues
            if f["location"]["start"]["value"] <= 90:
                has_tmh_n_terminus = True
                break  # Stop after finding the first one

    return (
        accession,
        organism_name,
        kingdom,
        protein_length,
        has_tmh_n_terminus,
        sequence
    )

def get_dataset(search_url, extract_function, output_file_name, header):
    """
    Runs the API search, processes all results, saves the formatted data to a TSV file,
    and returns a dictionary of accessions to sequences.
    """
    processed_entries_tsv = []
    processed_sequences = {}
    n_total = 0

    # Run the API call in batches
    for batch, total in get_batch(search_url):
        batch_json = json.loads(batch.text)
        for entry in batch_json["results"]:
            fields = extract_function(entry)
            if fields is not None:
                # Unpack all fields, the last one is the sequence
                *tsv_fields, sequence = fields
                accession = tsv_fields[0]
                
                processed_entries_tsv.append(tsv_fields)
                processed_sequences[accession] = sequence

            n_total += 1
        print(f"Processed {n_total} of {total} entries...")

    print(f"\nFinished processing. Found {len(processed_sequences)} proteins matching the criteria.")

    # Write the extracted metadata to a TSV file
    with open(output_file_name, "w") as ofs:
        # Write the new header row as specified
        ofs.write(header + "\n")
        for fields in processed_entries_tsv:
            print(*fields, sep="\t", file=ofs)

    print(f"Successfully saved data to {output_file_name}")
    return processed_sequences

def save_fasta_file(sequences_data, output_file_name):
    """
    Saves a dictionary of protein sequences to a FASTA file.
    """
    try:
        with open(output_file_name, "w") as f:
            for accession, sequence in sequences_data.items():
                f.write(f">{accession}\n")
                # Wrap sequence to 60 characters per line for standard FASTA format
                for i in range(0, len(sequence), 60):
                    f.write(sequence[i:i+60] + "\n")
        print(f"Successfully saved sequences to {output_file_name}")
    except IOError as e:
        print(f"Error writing FASTA file: {e}")

if __name__ == "__main__":
    # --- Positive Dataset ---
    print("--- Generating Positive Dataset ---")
    # URL updated to explicitly request the sequence and other necessary fields
    positive_url = "https://rest.uniprot.org/uniprotkb/search?format=json&query=%28%28fragment%3Afalse%29+AND+%28taxonomy_id%3A2759%29+AND+%28length%3A%5B40+TO+*%5D%29+AND+%28reviewed%3Atrue%29+AND+%28existence%3A1%29+AND+%28ft_signal_exp%3A*%29%29&size=500"
    positive_output_tsv = "positive_dataset_sp_cleavage.tsv"
    positive_output_fasta = "positive_dataset_sp_cleavage.fasta"
    positive_header = "Accession\tOrganism\tKingdom\tProtein_Length\tCleavage_Site"

    positive_sequences = get_dataset(positive_url, extract_fields_positive, positive_output_tsv, positive_header)
    save_fasta_file(positive_sequences, positive_output_fasta)

    # Read the TSV file into a pandas DataFrame and display a summary
    column_names_positive = [
        "Accession", "Organism", "Kingdom", "Protein_Length", "Cleavage_Site"
    ]
    df_positive = pd.read_csv(positive_output_tsv, sep='\t', names=column_names_positive, header=0)

    print(f"\nShape of the final positive dataset dataframe: {df_positive.shape}")
    print("\nFirst 10 rows of the positive dataset:")
    print(df_positive.head(10))
    print("\nKingdom distribution in the positive dataset:")
    print(df_positive['Kingdom'].value_counts())

    # --- Negative Dataset ---
    print("\n--- Generating Negative Dataset ---")
    # URL updated to explicitly request the sequence and other necessary fields
    negative_url = "https://rest.uniprot.org/uniprotkb/search?format=json&query=%28%28fragment%3Afalse%29+AND+%28reviewed%3Atrue%29+AND+%28existence%3A1%29+AND+%28taxonomy_id%3A2759%29+AND+%28length%3A%5B40+TO+*%5D%29+AND+%28%28cc_scl_term_exp%3ASL-0091%29+OR+%28cc_scl_term_exp%3ASL-0191%29+OR+%28cc_scl_term_exp%3ASL-0173%29+OR+%28cc_scl_term_exp%3ASL-0209%29+OR+%28cc_scl_term_exp%3ASL-0204%29+OR+%28cc_scl_term_exp%3ASL-0039%29%29+NOT+%28ft_signal%3A*%29%29&size=500"
    negative_output_tsv = "negative_dataset.tsv"
    negative_output_fasta = "negative_dataset.fasta"
    negative_header = "Accession\tOrganism\tKingdom\tProtein_Length\tTransmembrane_Helix_N_Terminus"

    negative_sequences = get_dataset(negative_url, extract_fields_negative, negative_output_tsv, negative_header)
    save_fasta_file(negative_sequences, negative_output_fasta)

    # Read the TSV file into a pandas DataFrame and display a summary
    column_names_negative = [
        "Accession", "Organism", "Kingdom", "Protein_Length", "Transmembrane_Helix_N_Terminus"
    ]
    df_negative = pd.read_csv(negative_output_tsv, sep='\t', names=column_names_negative, header=0)

    print(f"\nShape of the final negative dataset dataframe: {df_negative.shape}")
    print("\nFirst 10 rows of the negative dataset:")
    print(df_negative.head(10))
    print("\nDistribution of proteins with and without N-terminal TMH in the negative dataset:")
    print(df_negative['Transmembrane_Helix_N_Terminus'].value_counts())