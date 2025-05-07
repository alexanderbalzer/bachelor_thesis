import os
import pandas as pd
import logging
import shutil
from Bio import SeqIO
from datetime import datetime
from matplotlib import pyplot as plt
import numpy as np
import seaborn as sns

def format_species_name(name: str) -> str:
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


name = "Saccharomyces_cerevisiae"  # Caenorhabditis_elegans  Saccharomyces_cerevisiae  human
data = []  # Initialize an empty list to store the data
with open("pipeline/output/output_20250507_110227/" + name + "/" + name + "_filtered_by_GO_cleavable_mts.fasta", "r") as file:
    headers = []
    second_as = []
    for protein in SeqIO.parse(file, "fasta"):
        headers.append(protein.id)
        second_as.append(str(protein.seq)[1])
    # Create a pandas DataFrame with headers as the index and second_as as the column
    df_headers = pd.DataFrame({"Second_AS": second_as}, index=headers)
    

with open("pipeline/output/output_20250507_110227/"+ name + "/mitofates_for_" + name + ".cgi", "r") as file:
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
        if float(probability_of_MTS) >= 0.99:
            second_as = df_headers.loc[protein_id, "Second_AS"] if protein_id in df_headers.index else None
            if second_as is not None:
                # Append the data to the list
                data.append({"Protein_ID": protein_id, "Second_AS": second_as, "start_of_MTS": start_of_MTS, "Probability_of_MTS": probability_of_MTS, "length_of_MTS": length_of_MTS})
# Convert the list of dictionaries into a pandas DataFrame
df = pd.DataFrame(data)
print(len(df))
# Count the occurrences of each "Second_AS" at each "start_of_MTS"
heatmap_data = df.groupby(["start_of_MTS", "Second_AS"]).size().unstack(fill_value=0)
# Convert start_of_MTS to integers for chronological sorting
heatmap_data.index = heatmap_data.index.astype(int)
heatmap_data = heatmap_data.sort_index()
# Sum up all values in the heatmap_data
total_sum = heatmap_data.values.sum()
print(total_sum)
# Plot the heatmap
plt.figure(figsize=(10, 8))
sns.heatmap(heatmap_data, annot=True, fmt="d", cmap="viridis")
if name == "human" or name == "human_with_isoforms":
    biological_name = "Human"
else:
    biological_name = format_species_name(name)
plt.title("Heatmap for " + biological_name, fontstyle="italic")
colorbar = plt.gca().collections[0].colorbar
colorbar.set_label("Absolute count")
plt.xlabel("Second amino acid")
plt.ylabel("start of MTS")
plt.tight_layout()
plt.savefig("pipeline/output/output_20250507_110227/" + name + "/" + name + "_heatmap_start_of_MTS.png", dpi=300)

# Count the occurrences of each "Second_AS" for each "length_of_MTS"
heatmap_data = df.groupby(["length_of_MTS", "Second_AS"]).size().unstack(fill_value=0)
# Convert length_of_MTS to integers for chronological sorting
heatmap_data.index = heatmap_data.index.astype(int)
heatmap_data = heatmap_data.sort_index()
# Sum up all values in the heatmap_data
total_sum = heatmap_data.values.sum()
print(total_sum)
# Plot the heatmap
plt.figure(figsize=(10, 8))
sns.heatmap(heatmap_data, annot=True, fmt="d", cmap="PuBuGn")
if name == "human" or name == "human_with_isoforms":
    biological_name = "Human"
else:
    biological_name = format_species_name(name)
plt.title("Heatmap for " + biological_name, fontstyle="italic")
colorbar = plt.gca().collections[0].colorbar
colorbar.set_label("Absolute count")
plt.xlabel("Second amino acid")
plt.ylabel("length of MTS")
plt.tight_layout()
plt.savefig("pipeline/output/output_20250507_110227/" + name + "/" + name + "_heatmap_length_of_MTS.png", dpi=300)

