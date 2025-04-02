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
        "threshold": [""] * len(headers),
        "MTS-cleavable?": [""] * len(headers),
        "Sequence": sequences
    })
    
    return df

def parse_mts_cleavable_annotations(mts_file):
    """
    Parse MTS-cleavable annotations from a file.
    """
    cleavable = {}
    with open(mts_file, "r") as file:
        file_without_header = file.readlines()[2:]  
        for line in file_without_header:
            fields = line.strip().split("\t")
            if len(fields) < 3:
                continue
            protein_id = fields[0]
            cleavable[protein_id] = fields[2]
    logging.info(f"Parsed {len(cleavable)} MTS-cleavable annotations.")
    return cleavable

def parse_probability_of_mts (mts_file):
    """
    Parse probability of MTS-cleavable annotations from a file.
    """
    probability = {}
    with open(mts_file, "r") as file:
        file_without_header = file.readlines()[2:]  
        for line in file_without_header:
            fields = line.strip().split("\t")
            if len(fields) < 3:
                continue
            protein_id = fields[0]
            probability[protein_id] = float(fields[1])
    logging.info(f"Parsed {len(probability)} probability of MTS-cleavable annotations.")
    return probability


def add_mts_cleavable_to_dataframe(df, probability_of_mts):
    """
    Add MTS-cleavable information to the DataFrame based on the protein IDs.
    """
    for index, row in df.iterrows():
        protein_id = row["Header"]  
        if probability_of_mts[protein_id] > 0.9:
            df.at[index, "MTS-cleavable?"] = "Yes"
        elif probability_of_mts[protein_id] < 0.9:
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


name = ["human"]
name = name[0]

# Cleavable status can be "Yes" or "No"
cleavable = "No"
threshold = 0.9

fasta_file = "output files/filtered_proteins_by_GO_for_" + str(name) + ".fasta"
proteome = fasta_to_dataframe(fasta_file)

mts_cleavable_file = "input files/mito_fates_for_" + str(name) + ".cgi"

check_file_exists(fasta_file)
check_file_exists(mts_cleavable_file)

mts_cleavable = parse_mts_cleavable_annotations(mts_cleavable_file)
probability_of_mts = parse_probability_of_mts(mts_cleavable_file)

proteome = add_mts_cleavable_to_dataframe(proteome, probability_of_mts)

filtered_proteins = filter_proteins_by_mts(proteome, cleavable)
output_dir = "output files"


if cleavable == "No":
    output_file = "output files/filtered_proteins_no_cleavable_mts_for_" + str(name) + "_2.fasta"
if cleavable == "Yes":
    output_file = "output files/filtered_proteins_cleavable_mts_for_" + str(name) + "_2.fasta"


with open(output_file, "w") as f:
    for header, sequence in filtered_proteins:
        f.write(f">{header}\n{sequence}\n")