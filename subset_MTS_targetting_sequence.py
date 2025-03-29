import pandas as pd
from Bio import SeqIO
import os
import logging



logging.basicConfig(level=logging.INFO)

def check_file_exists(file_path):
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"File not found: {file_path}")


def fasta_to_dataframe(fasta_file):
    headers = []
    sequences = []

    for record in SeqIO.parse(fasta_file, "fasta"):
        headers.append(record.id)
        sequences.append(str(record.seq))
    
    df = pd.DataFrame({
        "Header": headers,
        "MTS-cleavable?": [""] * len(headers),    # Placeholder for GO terms
        "Sequence": sequences
    })
    
    return df

def parse_mts_cleavable_annotations(mts_file):
    """
    Parse MTS-cleavable annotations from a file.
    This function should be modified to read the actual MTS-cleavable data.
    """
    cleavable = {}
    with open(mts_file, "r") as file:
        file_without_header = file.readlines()[2:]  # Skip the header line
        for line in file_without_header:
            fields = line.strip().split("\t")
            if len(fields) < 3:
                continue
            protein_id = fields[0]
            cleavable[protein_id] = fields[2]
    logging.info(f"Parsed {len(cleavable)} MTS-cleavable annotations.")
    return cleavable


def add_mts_cleavable_to_dataframe(df, mts_cleavable):
    """
    Add MTS-cleavable information to the DataFrame based on the protein IDs.
    """
    for index, row in df.iterrows():
        protein_id = row["Header"]  # Extract protein ID from header
        mts_status = mts_cleavable.get(protein_id, "Unknown")  # Safely get the MTS status
        if mts_cleavable[protein_id] == "Possessing mitochondrial presequence":
            df.at[index, "MTS-cleavable?"] = "Yes"
        elif mts_cleavable[protein_id] == "No mitochondrial presequence":
            df.at[index, "MTS-cleavable?"] = "No"
        else:
            df.at[index, "MTS-cleavable?"] = "Unknown"
    return df

def filter_proteins_by_mts(dataframe, cleavable):
    cleavable_MTS = []
    for index, row in dataframe.iterrows():
        if row["MTS-cleavable?"] == cleavable:
            cleavable_MTS.append((row["Header"], row["Sequence"]))
    return cleavable_MTS



# Load the FASTA file
fasta_file = "output files/filtered_proteins.fasta"
proteome = fasta_to_dataframe(fasta_file)

# Load the MTS-cleavable information
mts_cleavable_file = "mts_cleavable_for_elegans.cgi"

# Check if files exist
check_file_exists(fasta_file)
check_file_exists(mts_cleavable_file)

# Parse the MTS-cleavable annotations
mts_cleavable = parse_mts_cleavable_annotations(mts_cleavable_file)

# Specify the MTS-cleavable status to filter by
# Options: "Yes" for cleavable, "No" for non-cleavable
cleavable = "Yes"

# Add MTS-cleavable information to the DataFrame
proteome = add_mts_cleavable_to_dataframe(proteome, mts_cleavable)

# Filter the proteins based on MTS-cleavable status
filtered_proteins = filter_proteins_by_mts(proteome, cleavable)
# Save the filtered proteins to a new FASTA file
output_dir = "output files"
os.makedirs(output_dir, exist_ok=True)  # Create the directory if it doesn't exist

output_file = os.path.join(output_dir, "filtered_proteins_mts_cleavable.fasta")
output_file = "output files/filtered_proteins_mts_cleavable.fasta"
with open(output_file, "w") as f:
    for header, sequence in filtered_proteins:
        f.write(f">{header}\n{sequence}\n")