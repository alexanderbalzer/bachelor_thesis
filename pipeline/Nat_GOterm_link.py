import pandas as pd
import numpy as np
import os
from Bio import SeqIO
import seaborn as sns
import matplotlib.pyplot as plt
from goatools.obo_parser import GODag
from goatools.base import download_go_basic_obo


# parse the fasta file and create a dataframe
def fasta_to_dataframe(fasta_file):
    data = []
    second_as = []
    protein_id = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        data.append(list(str(record.seq)))
        protein_id.append(record.id.split("|")[1])
    data = list(data)
    second_as = [seq[1] for seq in data]
    df = pd.DataFrame(second_as, columns=["2nd_amino_acid"])
    df["protein_id"] = protein_id
    os.makedirs("pipeline/protein_ids", exist_ok=True)
    with open("pipeline/protein_ids/protein_ids2.txt", "w") as f:
        for id in protein_id:
            f.write(f"{id}\n")
    return df

def parse_go_annotations(annotation_file, filter_by, filter):
    """
    Parse a GO annotation file to create a dictionary
    mapping protein IDs and their associated GO terms.
    """
    go_annotation = {}
    with open(annotation_file, "r") as file:
        for line in file:
            if line.startswith("!"):  
                continue
            fields = line.strip().split("\t")
            protein_id = fields[1]  #protein ID (column 2 in GAF)
            function = fields[8]  #function type (column 9 in GAF)
            if filter:
                if function != filter_by:  #filter for information
                    continue
            go_term = fields[4]     #GO term (column 5 in GAF)
            if protein_id not in go_annotation:
                go_annotation[protein_id] = []
            go_annotation[protein_id].append(go_term)
    return go_annotation

def add_go_terms_to_dataframe(df, go_annotation):
    """
    Add GO terms to the DataFrame based on the protein IDs.
    """
    for index, row in df.iterrows():
        protein_id = row["protein_id"]
        if protein_id in go_annotation:
            go_terms = go_annotation[protein_id]
            df.at[index, "GO_Term"] = ", ".join(go_terms)
        else:
            df.at[index, "GO_Term"] = "No GO term found"
    return df

def cluster_by_second_as(df):
    """
    Cluster the DataFrame by the second amino acid.
    """
    clusters = df.groupby("2nd_amino_acid").agg(list)
    clusters = clusters
    return clusters

def classify_natc_substrate(second_as):
    if second_as in ["A", "C", "T", "S", "V", "G"]:
        return "NatA/D"
    elif second_as in ["D", "E", "N", "Q"]:
        return "NatB"
    elif second_as in ["L", "I", "F", "Y", "K"]:
        return "NatC/E"
    else:
        return "Other"
    
def get_go_aspect(go_id):
    """
    Get the aspect of a GO term (biological process, molecular function, or cellular component).
    """
    if go_id in go_dag:
        go_term = go_dag[go_id]
        return go_term.name
    else:
        return "Unknown"

def format_species_name(name: str) -> str:
    # Teile den Namen anhand des Unterstrichs
    parts = name.split("_")
    if name == "human":
        return "Human"
    if name == "human_with_isoforms":
        return "Human with isoforms"
    if len(parts) != 2:
        raise ValueError("Name muss genau ein Unterstrich enthalten (Gattung_Art)")
    
    genus, species = parts
    # Kürze den Gattungsnamen auf den ersten Buchstaben + Punkt
    short_genus = genus[0] + "."
    
    # Setze alles in kursiv (z. B. für Markdown oder HTML)
    formatted = f"{short_genus} {species}"
    return formatted


    
obo_path = "pipeline/go.obo"
if not os.path.exists(obo_path):
    # Download the latest OBO file
    obo_path = download_go_basic_obo()  # Downloads the latest OBO file

# Load the GO DAG
go_dag = GODag(obo_path)



filter = True
threshold = 5  # Minimum number of occurrences for a GO term to be included in the heatmap
filter_by = "C"  
# filter_by = "C"  # Filter for cellular component
# filter_by = "F"  # Filter for function
# filter_by = "P"  # Filter for process
name = "human" # Saccharomyces_cerevisiae Caenorhabditis_elegans human



def run(organism_names, input_dir, ouput_dir):
    for i, name in enumerate(organism_names, start=1):
        print(f"Processing organism: {name}")
        output_dir_per_organism = ouput_dir + "/" + name
        fasta = output_dir_per_organism + "/" + name + "_filtered_by_GO_cleavable_mts.fasta"
        goa = os.path.join(input_dir, f"{name}.goa")
        df = fasta_to_dataframe(fasta)
        go_annotation = parse_go_annotations(goa, filter_by, filter)
        df = add_go_terms_to_dataframe(df, go_annotation)
        df = df.drop(columns=["protein_id"])
        # Create a new DataFrame with one-hot encoding for GO terms
        go_term_counts = df["GO_Term"].str.get_dummies(sep=", ")

        # Add a column "NatC_Substrate" based on the second amino acid
        df["NatC_Substrate"] = df["2nd_amino_acid"].apply(classify_natc_substrate)

        # Add the counts to the original DataFrame
        df = pd.concat([df, go_term_counts], axis=1)

        # Group by the second amino acid and sum the counts for each GO term
        go_term_summary = df.groupby("NatC_Substrate")[go_term_counts.columns].sum()
        # Remove columns with less than 4 occurrences for any amino acid
        go_term_summary = go_term_summary.loc[:, (go_term_summary >= threshold).any()]
        # transpose the DataFrame for better visualization
        go_term_summary = go_term_summary.transpose()
            # Replace the GO term "GO:0005739" with its sum and move it to the bottom
        if "GO:0005739" in go_term_summary.index:
            go_term_summary.loc["sum"] = go_term_summary.loc["GO:0005739"]
            go_term_summary = go_term_summary.drop("GO:0005739")
        print(go_term_summary)
        
        if filter_by == "C":
            function = "cellular component"
        elif filter_by == "F":
            function = "function"
        elif filter_by == "P":
            function = "process"
        else:
            function = "all_functions"
        # exchange GO terms with their aspects
        go_term_summary.index = [get_go_aspect(go_id) for go_id in go_term_summary.index]
        if "Unknown" in go_term_summary.index:
            go_term_summary = go_term_summary.drop("Unknown")

        # Create a heatmap
        plt.figure(figsize=(12, 8))
        if go_term_summary.empty:
            print(f"No data available for {name} with filter {filter_by}.")
            continue
        sns.heatmap(go_term_summary, cmap="YlGnBu", annot=True, fmt="d")
        if name == "human":
            biological_name = "Human"
        else:
            biological_name = format_species_name(name)
        plt.title(biological_name, fontstyle="italic")
        colorbar = plt.gca().collections[0].colorbar
        colorbar.set_label("Absolute count")
        plt.xlabel("Nat substrate")
        plt.ylabel(function)
        plt.xticks(rotation=90)
        plt.tight_layout()
        os.makedirs(output_dir_per_organism, exist_ok=True)
        plt.savefig(os.path.join(output_dir_per_organism, f"heatmap_{name}_{function}.png"))
        plt.close()
    
if __name__ == "__main__":
    # Define the names of the organisms
    organism_names = [
        "Geotrichum_candidum", "Drosophila_Melanogaster", "Arabidopsis_thaliana", 
        "Lachancea_thermotolerans", "human_with_isoforms", "human", "Clavispora_lusitaniae", 
        "Mus_musculus", "Caenorhabditis_elegans", "Candida_glabrata", "Schizosaccharomyces_pombe", 
        "Debaryomyces_hansenii", "Yarrowia_lipolytica", "Saccharomyces_cerevisiae", 
        "Zygosaccharomyces_rouxii", "Physcomitrium_patens", "Scheffersomyces_stipitis"]
    output_dir = "pipeline/output/output_20250506_171728"
    input_dir = "pipeline/input"
    run(organism_names, input_dir, output_dir)
