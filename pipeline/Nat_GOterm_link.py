import pandas as pd
import numpy as np
import os
from Bio import SeqIO
import seaborn as sns
import matplotlib.pyplot as plt
from goatools.obo_parser import GODag
from goatools.base import download_go_basic_obo
from time import sleep


def fasta_to_dataframe(fasta_file):
    """
    Parse a FASTA file and extract the second amino acid and protein ID.
    """
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
    
def get_go_aspect(go_id, go_dag):
    """
    Get the aspect of a GO term.
    """
    if go_id in go_dag:
        go_term = go_dag[go_id]
        return go_term.name
    else:
        return "Unknown"

def format_species_name(name: str) -> str:
    # Teile den Namen anhand des Unterstrichs
    parts = name.split("_")
    if name == "Homo_sapiens_isoforms":
        return "H. Sapiens with isoforms"
    if len(parts) != 2:
        raise ValueError("Name has to contain two words.")
    genus, species = parts
    short_genus = genus[0] + "."
    formatted = f"{short_genus} {species}"
    return formatted


def load_obo():
    '''
    load the latest obo file to exchange the GO terms with their aspects
    '''
    obo_path = "pipeline/go.obo"
    if not os.path.exists(obo_path):
        obo_path = download_go_basic_obo()  
    # Load the GO DAG
    go_dag = GODag(obo_path)
    return go_dag



def run(organism_names, input_dir, ouput_dir):
    """
    # plot the frequency of the Nat substrates for each GO term inside the mitochondria for each organism
    # and save the heatmap as a pdf file
    Args:
        organism_names (list): List of organism names.
        input_dir (str): Directory containing the input files.
        ouput_dir (str): Directory to save the output files.
    """
    filter = True
    absolute_threshold = 5
    filter_by = "C" 
    nat_or_aa = "aa" # nat or aa
    # Load the GO DAG
    go_dag = load_obo()
    for i, name in enumerate(organism_names, start=1):
        print(f"Creating Nat - GO term heatmap for organism: {name}")

        # specify the input files for cleavable MTS
        output_dir_per_organism = ouput_dir + "/" + name
        fasta = output_dir_per_organism + "/" + name + "_filtered_by_GO_cleavable_mts.fasta"
        goa = os.path.join(input_dir, f"{name}.goa")

        # parse the fasta file and the goa file and create a dataframe
        df = fasta_to_dataframe(fasta)
        go_annotation = parse_go_annotations(goa, filter_by, filter)
        df = add_go_terms_to_dataframe(df, go_annotation)

        # exchange the go terms with their aspects
        df["GO_Term"] = df["GO_Term"].apply(lambda x: ", ".join([get_go_aspect(go_id, go_dag) for go_id in x.split(", ")]))
        # delte duplicates in the GO_Term column
        df["GO_Term"] = df["GO_Term"].apply(lambda x: ", ".join(sorted(set(x.split(", ")))))
        # delete the term mitochondrion from the GO_Term column
        df["GO_Term"] = df["GO_Term"].apply(lambda x: ", ".join([term for term in x.split(", ") if "mitochondrion" not in term]))
        # delete all terms that don't contain "mitochondr"
        df["GO_Term"] = df["GO_Term"].apply(lambda x: ", ".join([term for term in x.split(", ") if "mitochondr" in term]))
        # delete all terms that are empty
        df = df[df["GO_Term"] != ""]
        print(len(df))
        # save the dataframe with the GO terms
        df.to_csv(os.path.join(output_dir_per_organism, f"{name}_go_terms.csv"), index=False)

        # Explode GO_Term so each row has one GO term, then group
        df_exploded = df.copy()
        df_exploded = df_exploded.assign(GO_Term=df_exploded["GO_Term"].str.split(", ")).explode("GO_Term")
        go_term_protein_ids = df_exploded.groupby("GO_Term")["protein_id"].apply(list).reset_index()
        go_term_protein_ids.columns = ["GO_Term", "Protein_IDs"]
        # add the number of proteins per go term
        go_term_protein_ids["Number_of_Proteins"] = go_term_protein_ids["Protein_IDs"].apply(lambda x: len(x))
        # Save the DataFrame to a CSV file
        go_term_protein_ids.to_csv(os.path.join(output_dir_per_organism, f"{name}_go_term_protein_ids.csv"), index=False)
 
        df = df.drop(columns=["protein_id"])
        # Create a new DataFrame with one-hot encoding for GO terms
        go_term_counts = df["GO_Term"].str.get_dummies(sep=", ")

        # Add a column "NatC_Substrate" based on the second amino acid
        df["NatC_Substrate"] = df["2nd_amino_acid"].apply(classify_natc_substrate)
        # Add the counts to the original DataFrame
        df = pd.concat([df, go_term_counts], axis=1)


        # Group by the second amino acid and sum the counts for each GO term
        if nat_or_aa == "nat":
            go_term_summary = df.groupby("NatC_Substrate")[go_term_counts.columns].sum()
        else:
            go_term_summary = df.groupby("2nd_amino_acid")[go_term_counts.columns].sum()
        

        # transpose the DataFrame for better visualization
        go_term_summary = go_term_summary.transpose()

        # Remove rows with less than 4 proteins
        for row in go_term_summary.index:
            if go_term_summary.loc[row].sum() < absolute_threshold:
                go_term_summary = go_term_summary.drop(row)

        # append the sum of each row to the index
        go_term_summary.index = go_term_summary.index + ", (n=" + go_term_summary.sum(axis=1).astype(str) + ")"


        # normalize the counts rowwise
        go_term_summary = go_term_summary.div(go_term_summary.sum(axis=1), axis=0)
        
        if filter_by == "C":
            function = "cellular component"
        elif filter_by == "F":
            function = "function"
        elif filter_by == "P":
            function = "process"
        else:
            function = "all_functions"
        if "Unknown" in go_term_summary.index:
            go_term_summary = go_term_summary.drop("Unknown")

        # Create a heatmap
        plt.figure(figsize=(12, 8))
        if go_term_summary.empty:
            print(f"No data available for {name} with filter {filter_by}.")
            continue

        # Reorder the columns to match the specified amino acid order if the x axis is amino acids
        if nat_or_aa == "aa":
            amino_acid_order = ["D", "E", "N", "Q", "Y", "H", "K", "R", "M", "L", "F", "I", "W", "S", "A", "T", "C", "P", "G", "V"]
            reordered_indices = [go_term_summary.columns.tolist().index(aa) for aa in amino_acid_order if aa in go_term_summary.columns]
            go_term_summary = go_term_summary.iloc[:, reordered_indices]
        # Create the heatmap
        sns.heatmap(go_term_summary, cmap="YlGnBu", annot=False, fmt=".2f")
        if name == "human":
            biological_name = "Human"
        else:
            biological_name = format_species_name(name)
        plt.title(biological_name, fontstyle="italic")
        colorbar = plt.gca().collections[0].colorbar
        colorbar.set_label("Relative abundance")
        if nat_or_aa == "nat":
            plt.xlabel("NatC Substrates")
        else:
            plt.xlabel("Amino acids")
        plt.ylabel(function)
        plt.xticks(rotation=0)
        plt.tight_layout()
        os.makedirs(output_dir_per_organism, exist_ok=True)
        plt.savefig(os.path.join(output_dir_per_organism, f"heatmap_{name}_{function}_noncleavable.pdf"))
        # save the matrix as a csv file
        go_term_summary.to_csv(os.path.join(output_dir_per_organism, f"heatmap_{name}_{function}_noncleavable.csv"))
        plt.close()
    
if __name__ == "__main__":
    # Define the names of the organisms
    organism_names = [
    "Homo_sapiens", "Homo_sapiens_isoforms", "Mus_musculus", "Dario_rerio", "Daphnia_magna", 
    "Caenorhabditis_elegans", "Drosophila_Melanogaster", "Arabidopsis_thaliana", 
    "Physcomitrium_patens", "Chlamydomonas_reinhardtii", 
    "Candida_glabrata", "Saccharomyces_cerevisiae", "Zygosaccharomyces_rouxii"]
    output_dir = "pipeline/output/output_20250513_105122"
    input_dir = "pipeline/input"
    filter = True
    absolute_threshold = 5
    relative_threshold = 0.5  # Minimum number of occurrences for a GO term to be included in the heatmap
    filter_by = "C" 
    run(organism_names, input_dir, output_dir)
