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
        sequence = row["MTS_Sequence"]
        new_sequence = "".join([property_mapping.get(aa, aa) for aa in sequence])
        df.at[index, "MTS_Sequence"] = new_sequence
    return df


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
            # only include proteins with a G as second amino acid
#            df_headers = df_headers[df_headers["Second_AS"] == "L"]
            
        if start == "MTS":
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

        if start == "MTS":
            # Convert the list of dictionaries into a pandas DataFrame
            df = pd.DataFrame(data)
            df["MTS_Sequence"] = df.apply(
            lambda row: row["Sequence"][int(row["start_of_MTS"]) - 1:int(row["start_of_MTS"]) - 1 + int(row["length_of_MTS"])],
            axis=1
            )
        elif start == "beginning":
            df = df_headers
            df["MTS_Sequence"] = df.apply(
            lambda row: row["Sequence"][1:1 + int(50)],
            axis=1
            )
        else:
            raise ValueError("Invalid start value. Use 'MTS' or 'beginning'.")

        df = replace_amino_acids_with_properties(df)
        # Generate the frequency matrix
        frequency_matrix = create_frequency_matrix(df["MTS_Sequence"])
        frequency_matrix = frequency_matrix.fillna(0)
        # Create a color scheme for the logoplot
        # Define a custom color scheme
        return frequency_matrix, df
custom_color_scheme2 = {
    'L': 'green', 'F': 'green', 'I': 'green', 'V': 'green', 'W': 'green', 'Y': 'green', 'M': 'green', 'C': 'green', 'A': 'green',
    'R': 'blue', 'K': 'blue', 'H': 'blue',
    'S': 'lightblue', 'T': 'lightblue', 'N': 'lightblue', 'Q': 'lightblue',
    'P': 'yellow', 'G': 'yellow',
    'D': 'gray', 'E': 'gray', 'X': 'gray',
    '-': 'white'
}
# Define a custom color scheme
custom_color_scheme = {
    'H': 'green',
    'B': 'blue',
    'P': 'lightblue',
    'X': 'yellow',
    '-': 'black'
}
# Create a legend for the color scheme
legend_elements = [
    Patch(facecolor='green', edgecolor='black', label='Hydrophobic'),
    Patch(facecolor='blue', edgecolor='black', label='Basic'),
    Patch(facecolor='lightblue', edgecolor='black', label='Polar'),
    Patch(facecolor='yellow', edgecolor='black', label='Secondary Structure Breaker')
#    Patch(facecolor='gray', edgecolor='black', label='Not Relevant')
]
def run_MTS_and_start(organism_names, output_dir):
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
        frequency_matrix_mts = frequency_matrix_mts.iloc[:9, :]  # Cut after the 10th row
        logo_mts = Logo(frequency_matrix_mts, color_scheme=custom_color_scheme, ax=ax_mts)  # Limit to first 8 positions
        ax_mts.set_title(f"Logoplot of MTS Sequences for {format_species_name(name)}")
        ax_mts.set_xticks(range(9))  # Set x-axis ticks for positions 1 to 8
        ax_mts.set_xticklabels(range(1, 10))  # Set x-axis labels for positions 1 to 8

        # Create the beginning sequence logoplot
        ax_beginning = axes[1]
        frequency_matrix_beginning = frequency_matrix_beginning.iloc[:9, :]  # Limit to first 9 positions
        logo_beginning = Logo(frequency_matrix_beginning, color_scheme=custom_color_scheme, ax=ax_beginning)  # Limit to first 9 positions
        ax_beginning.set_title(f"Logoplot of Beginning Sequences for {format_species_name(name)}")
        ax_beginning.set_xticks(range(9))  # Set x-axis ticks for positions 1 to 8
        ax_beginning.set_xticklabels(range(1, 10))  # Set x-axis labels for positions 1 to 8

        # Add a legend to the first subplot
        ax_mts.legend(handles=legend_elements, title="Amino Acid Properties", loc='upper right')

        # Adjust layout and save the figure
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir_per_organism, f"logoplot_{name}_L_both.png"))
        plt.close(fig)
        # save whole data to fasta file
        with open(os.path.join(output_dir_per_organism, f"MTS_sequences_{name}.fasta"), "w") as fasta_file:
            for index, row in data.iterrows():
                fasta_file.write(f">{row['Protein_ID']}\n+{row['MTS_Sequence']}\n<{row['Sequence']}\n*{row['start_of_MTS']}\n")


def run_start(organism_names, output_dir):
    for name in organism_names:
        print(f"Processing organism: {name}")
        output_dir_per_organism = output_dir + "/" + name 
        # Generate frequency matrix for beginning sequences
        frequency_matrix_beginning, unused = generate_frequency_matrix(name, "beginning", output_dir_per_organism)
        # Create a figure for the beginning sequence logoplot
        fig, ax = plt.subplots(figsize=(15, 5))
        frequency_matrix_beginning = frequency_matrix_beginning.iloc[:9, :]
        logo_beginning = Logo(frequency_matrix_beginning, color_scheme=custom_color_scheme, ax=ax)
        ax.set_title(f"Logoplot of Beginning Sequences for {format_species_name(name)}")
        ax.set_xticks(range(9))  # Set x-axis ticks for positions 1 to 8
        ax.set_xticklabels(range(1, 10))
        ax.legend(handles=legend_elements, title="Amino Acid Properties", loc='upper right')
        # add y axis grid lines
        ax.yaxis.grid(True, linestyle='--', alpha=0.7)

        # Adjust layout and save the figure
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir_per_organism, f"logoplot_beginning_{name}.png"))
        plt.close(fig)


if __name__ == "__main__":
    # Define the names of the organisms
    organism_names = [
        "Geotrichum_candidum", "Drosophila_Melanogaster", "Arabidopsis_thaliana", 
        "Lachancea_thermotolerans", "human_with_isoforms", "human", "Clavispora_lusitaniae", 
        "Mus_musculus", "Caenorhabditis_elegans", "Candida_glabrata", "Schizosaccharomyces_pombe", 
        "Debaryomyces_hansenii", "Yarrowia_lipolytica", "Saccharomyces_cerevisiae", 
        "Zygosaccharomyces_rouxii", "Physcomitrium_patens", "Scheffersomyces_stipitis"]
    output_dir = "pipeline/output/output_20250509_152503"

    run_start(organism_names, output_dir)



