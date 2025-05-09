import Bio
import pandas as pd
import numpy as np
import os
from Bio import SeqIO
import seaborn as sns
import matplotlib.pyplot as plt
from goatools.obo_parser import GODag   

def replace_amino_acids_with_properties(df):
    """
    replace the sequences in the DataFrame based on their properties.
    """
    # Define the mapping of amino acids to properties
    property_mapping = {
        'L': 'H', 'F': 'H', 'I': 'H', 'V': 'H', 'W': 'H', 
        'Y': 'H', 'M': 'H', 'C': 'H', 'A': 'H',
        'R': 'B', 'K': 'B', 'H': 'B',
        'S': 'P', 'T': 'P', 'N': 'P', 'Q': 'P',
        'P': 'X', 'G': 'X', 'E': 'X', 'D': 'X',
        '-': 'X'
    }
    # Replace the sequences in the DataFrame based on their properties
    for index, row in df.iterrows():
        sequence = row["sequence"]
        new_sequence = "".join([property_mapping.get(aa, aa) for aa in sequence])
        # add the new sequence to the DataFrame if the first 5 positions are H
        if new_sequence[:5] == "HHHHH":
            df.at[index, "sequence"] = new_sequence
    return df

def fasta_to_dataframe(fasta_file):
    data = []
    protein_id = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        data.append(str(record.seq))
        protein_id.append(record.id)
    data = list(data)
    df = pd.DataFrame(data, columns=["sequence"])
    # replace the sequences in the DataFrame based on their properties
    df = replace_amino_acids_with_properties(df)
    df["protein_id"] = protein_id
    return df

def run(organism_list, output_dir):
    for i, name in enumerate(organism_list, start=1):
        output_dir_per_organism = output_dir + "/" + name
        print(f"Processing organism: {name}")
        with open(output_dir_per_organism + "/" + name + "_filtered_by_GO_cleavable_mts.fasta", "r") as file:
            df = fasta_to_dataframe(file)
        # write the modified DataFrame to a new file
        output_file = os.path.join(output_dir_per_organism, f"{name}_modified.fasta")
        with open(output_file, "w") as f:
            for index, row in df.iterrows():
                f.write(f">{row['protein_id']}\n")
                f.write(f"{row['sequence']}\n")
        print(f"Modified sequences saved to {output_file}")
    



            

if __name__ == "__main__":
    output_dir = "pipeline/output/output_20250508_153218_ER"
    # Define the names of the organisms
    organism_list = [
        "Geotrichum_candidum", "Drosophila_Melanogaster", "Arabidopsis_thaliana", 
        "Lachancea_thermotolerans", "human_with_isoforms", "human", "Clavispora_lusitaniae", 
        "Mus_musculus", "Caenorhabditis_elegans", "Candida_glabrata", "Schizosaccharomyces_pombe", 
        "Debaryomyces_hansenii", "Yarrowia_lipolytica", "Saccharomyces_cerevisiae", 
        "Zygosaccharomyces_rouxii", "Physcomitrium_patens", "Scheffersomyces_stipitis"]
    run(organism_list, output_dir)