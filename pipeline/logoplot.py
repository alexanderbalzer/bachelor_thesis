import os
import pandas as pd
import logging
import shutil
from Bio import SeqIO
from datetime import datetime
from matplotlib import pyplot as plt
import numpy as np
import seaborn as sns
from logomaker import Logo
from collections import Counter
from matplotlib.patches import Patch
from itertools import islice

def format_species_name(name: str) -> str:
    if name == "human":
        return "human"
    elif name == "human_with_isoforms":
        return "human with isoforms"
    # Teile den Namen anhand des Unterstrichs
    parts = name.split("_")
    if len(parts) != 2:
        raise ValueError("Name muss genau ein Unterstrich enthalten (Gattung_Art)")
    
    genus, species = parts
    # Kürze den Gattungsnamen auf den ersten Buchstaben + Punkt
    short_genus = genus[0] + "."
    
    # Setze alles in kursiv (z. B. für Markdown oder HTML)
    formatted = f"{short_genus} {species}"
    return formatted

# Create a frequency matrix for the MTS sequences
def create_frequency_matrix(sequences):
    max_length = max(len(seq) for seq in sequences)
    frequency_matrix = []
    for i in range(max_length):
        column = [seq[i] if i < len(seq) else '-' for seq in sequences]
        counts = Counter(column)
        total = sum(counts.values())
        frequency_matrix.append({key: counts[key] / total for key in counts})
    return pd.DataFrame(frequency_matrix).fillna(0)

def generate_frequency_matrix(name, start, output_dir_per_organism):
        data = []  # Initialize an empty list to store the data
        with open(output_dir_per_organism + "/" + name + "_filtered_by_GO_cleavable_mts.fasta", "r") as file:
            headers = []
            second_as = []
            sequences = []
            for protein in SeqIO.parse(file, "fasta"):
                headers.append(protein.id)
                sequences.append(str(protein.seq))
                second_as.append(str(protein.seq)[1])
            # Create a pandas DataFrame with headers as the index and second_as as the column
            df_headers = pd.DataFrame({"Second_AS": second_as, "Sequence": sequences}, index=headers)
            

        with open(output_dir_per_organism + "/" + "mitofates_for_" + name + ".cgi", "r") as file:
            lines = file.readlines()
            for i, line in enumerate(lines):
                if line.startswith("!") or i == 0:  # Skip header line
                    continue
                fields = line.strip().split("\t")
                protein_id = fields[0]
                probability_of_MTS = fields[1]
                position_of_MTS = fields[6]
                start_of_MTS = position_of_MTS.strip().split("-")[0]
                length_of_MTS = float(position_of_MTS.strip().split("-")[1]) - float(position_of_MTS.strip().split("-")[0])
                
                # Append the data for this line to the list if the probability is above the threshold
                if float(probability_of_MTS) >= 0.9:
                    second_as = df_headers.loc[protein_id, "Second_AS"] if protein_id in df_headers.index else None
                    sequence = df_headers.loc[protein_id, "Sequence"] if protein_id in df_headers.index else None
                    if second_as is not None:
                        # Append the data to the list
                        data.append({"Protein_ID": protein_id, "Sequence": sequence, "Second_AS": second_as, "start_of_MTS": start_of_MTS, "Probability_of_MTS": probability_of_MTS, "length_of_MTS": length_of_MTS})
        # Convert the list of dictionaries into a pandas DataFrame
        df = pd.DataFrame(data)
        # Create a new column for the desired sequence
        if start == "MTS":
            df["MTS_Sequence"] = df.apply(
            lambda row: row["Sequence"][int(row["start_of_MTS"]) - 1:int(row["start_of_MTS"]) - 1 + int(row["length_of_MTS"])],
            axis=1
            )
        elif start == "beginning":
            df["MTS_Sequence"] = df.apply(
            lambda row: row["Sequence"][1:1 + int(row["length_of_MTS"])],
            axis=1
            )
        else:
            raise ValueError("Invalid start value. Use 'MTS' or 'beginning'.")


        # Generate the frequency matrix
        frequency_matrix = create_frequency_matrix(df["MTS_Sequence"])
        frequency_matrix = frequency_matrix.fillna(0)
        # Create a color scheme for the logoplot
        # Define a custom color scheme
        return frequency_matrix, df
custom_color_scheme = {
    'L': 'green', 'F': 'green', 'I': 'green', 'V': 'green', 'W': 'green', 'Y': 'green', 'M': 'green', 'C': 'green', 'A': 'green',
    'R': 'blue', 'K': 'blue', 'H': 'blue',
    'S': 'lightblue', 'T': 'lightblue', 'N': 'lightblue', 'Q': 'lightblue',
    'P': 'yellow', 'G': 'yellow',
    'D': 'gray', 'E': 'gray', 'X': 'gray',
    '-': 'white'

}
# Create a legend for the color scheme
legend_elements = [
    Patch(facecolor='green', edgecolor='black', label='Hydrophobic'),
    Patch(facecolor='blue', edgecolor='black', label='Basic'),
    Patch(facecolor='lightblue', edgecolor='black', label='Polar'),
    Patch(facecolor='yellow', edgecolor='black', label='Secondary Structure Breaker')
#    Patch(facecolor='gray', edgecolor='black', label='Not Relevant')
]
def run(organism_names, output_dir):
    for name in organism_names:
        print(f"Processing organism: {name}")
        output_dir_per_organism = output_dir + "/" + name 
        # Generate frequency matrices for both MTS and beginning sequences
        frequency_matrix_mts, data = generate_frequency_matrix(name, "MTS", output_dir_per_organism)
        frequency_matrix_beginning, unused = generate_frequency_matrix(name, "beginning", output_dir_per_organism)

        # Create a figure with two subplots (one for MTS, one for beginning sequence)
        fig, axes = plt.subplots(2, 1, figsize=(15, 10))  # Two rows, one column

        # Create the MTS logoplot
        ax_mts = axes[0]
        logo_mts = Logo(frequency_matrix_mts, color_scheme=custom_color_scheme, ax=ax_mts)
        ax_mts.set_title(f"Logoplot of MTS Sequences for {format_species_name(name)}")

        # Create the beginning sequence logoplot
        ax_beginning = axes[1]
        logo_beginning = Logo(frequency_matrix_beginning, color_scheme=custom_color_scheme, ax=ax_beginning)
        ax_beginning.set_title(f"Logoplot of Beginning Sequences for {format_species_name(name)}")

        # Add a legend to the first subplot
        ax_mts.legend(handles=legend_elements, title="Amino Acid Properties", loc='upper right')

        # Adjust layout and save the figure
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir_per_organism, f"logoplot_{name}.png"))
        plt.close(fig)
        # save whole data to fasta file
        with open(os.path.join(output_dir_per_organism, f"MTS_sequences_{name}.fasta"), "w") as fasta_file:
            for index, row in data.iterrows():
                fasta_file.write(f">{row['Protein_ID']}\n+{row['MTS_Sequence']}\n<{row['Sequence']}\n*{row['start_of_MTS']}\n")


if __name__ == "__main__":
    # Define the names of the organisms
    organism_names = [
        "Geotrichum_candidum", "Drosophila_Melanogaster", "Arabidopsis_thaliana", 
        "Lachancea_thermotolerans", "human_with_isoforms", "human", "Clavispora_lusitaniae", 
        "Mus_musculus", "Caenorhabditis_elegans", "Candida_glabrata", "Schizosaccharomyces_pombe", 
        "Debaryomyces_hansenii", "Yarrowia_lipolytica", "Saccharomyces_cerevisiae", 
        "Zygosaccharomyces_rouxii", "Physcomitrium_patens", "Scheffersomyces_stipitis"]
    output_dir = "pipeline/output/output_20250430_143623"

    run(organism_names, output_dir)
    


