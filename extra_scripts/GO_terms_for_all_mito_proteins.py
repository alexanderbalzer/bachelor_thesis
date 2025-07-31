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
    ''' Parse a FASTA file and create a DataFrame with protein IDs and their second amino acids. '''
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
    if name == "Homo_sapiens_isoforms":
        return "H. Sapiens with isoforms"
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
threshold = 0.5  # Minimum number of occurrences for a GO term to be included in the heatmap
filter_by = "C"  
# filter_by = "C"  # Filter for cellular component
# filter_by = "F"  # Filter for function
# filter_by = "P"  # Filter for process
name = "human" # Saccharomyces_cerevisiae Caenorhabditis_elegans human



def run(organism_names, input_dir, ouput_dir, relative_threshold, filter, filter_by, absolute_threshold):
    """
    Fetches all GO terms for the specified organisms, processes them, and generates heatmaps 
    showing the relative usage of NATC substrates across different GO terms divided into MTS bearing and non-MTS bearing proteins.
    Input:
    - organism_names: List of organism names to process.
    - input_dir: Directory containing the input files.
    - ouput_dir: Directory to save the output files.
    - relative_threshold: Minimum number of occurrences for a GO term to be included in the heatmap.
    - filter: Boolean indicating whether to filter the GO terms.
    - filter_by: The type of GO term to filter by (C, F, P).
    - absolute_threshold: Minimum number of occurrences for a GO term to be included in the heatmap.
    """
    for i, name in enumerate(organism_names, start=1):
        arrays = []
        for z in range(1, 3):
            print(f"Processing organism: {name}")
            output_dir_per_organism = ouput_dir + "/" + name
            if z == 1:
                fasta = output_dir_per_organism + "/" + name + "_filtered_by_GO_cleavable_mts.fasta"
            elif z == 2:
                fasta = output_dir_per_organism + "/" + name + "_filtered_proteins_by_GO_noncleavable_mts.fasta"
            goa = os.path.join(input_dir, f"{name}.goa")
            df = fasta_to_dataframe(fasta)
            go_annotation = parse_go_annotations(goa, filter_by, filter)
            df = add_go_terms_to_dataframe(df, go_annotation)

            # exchange the go terms with their aspects
            df2 = df.copy()
            df2["GO_Term"] = df2["GO_Term"].apply(lambda x: ", ".join([get_go_aspect(go_id) for go_id in x.split(", ")]))
            # delte duplicates in the GO_Term column
            df2["GO_Term"] = df2["GO_Term"].apply(lambda x: ", ".join(sorted(set(x.split(", ")))))
            # delete the term mitochondrion from the GO_Term column
            #df2["GO_Term"] = df2["GO_Term"].apply(lambda x: ", ".join([term for term in x.split(", ") if "mitochondrion" not in term]))
            # delete all terms that don't contein "mitochondr"
            df2["GO_Term"] = df2["GO_Term"].apply(lambda x: ", ".join([term for term in x.split(", ") if "mitochondr" in term]))
            # save the dataframe with the GO terms
            df2.to_csv(os.path.join(output_dir_per_organism, f"{name}_go_terms.csv"), index=False)
            # delete all terms that are empty
            df2 = df2[df2["GO_Term"] != ""]
            print(len(df2))

            df = df2
            df = df.drop(columns=["protein_id"])
            # Create a new DataFrame with one-hot encoding for GO terms
            go_term_counts = df["GO_Term"].str.get_dummies(sep=", ")

            # Add a column "NatC_Substrate" based on the second amino acid
            df["NatC_Substrate"] = df["2nd_amino_acid"].apply(classify_natc_substrate)
            # Add the counts to the original DataFrame
            df = pd.concat([df, go_term_counts], axis=1)


            # Group by the second amino acid and sum the counts for each GO term
            go_term_summary = df.groupby("NatC_Substrate")[go_term_counts.columns].sum()

            # transpose the DataFrame for better visualization
            go_term_summary = go_term_summary.transpose()


            # exchange GO terms with their aspects
            #go_term_summary.index = [get_go_aspect(go_id) for go_id in go_term_summary.index]


            # Remove rows with less than 4 proteins
            for row in go_term_summary.index:
                if go_term_summary.loc[row].sum() < absolute_threshold:
                    go_term_summary = go_term_summary.drop(row)

            # Remove columns with less than 4 occurrences for any amino acid
            #go_term_summary = go_term_summary.loc[:, (go_term_summary >= relative_threshold).any()]
            
            if filter_by == "C":
                function = "cellular component"
            elif filter_by == "F":
                function = "function"
            elif filter_by == "P":
                function = "process"
            else:
                function = "all_functions"
            '''if "Unknown" in go_term_summary.index:
                go_term_summary = go_term_summary.drop("Unknown")'''

            # Create a heatmap
            plt.figure(figsize=(12, 8))
            if go_term_summary.empty:
                print(f"No data available for {name} with filter {filter_by}.")
                continue

            # Reorder the columns to match the specified amino acid order
            '''amino_acid_order = ["D", "E", "N", "Q", "Y", "H", "K", "R", "M", "L", "F", "I", "W", "S", "A", "T", "C", "P", "G", "V"]
            reordered_indices = [go_term_summary.columns.tolist().index(aa) for aa in amino_acid_order if aa in go_term_summary.columns]
            go_term_summary = go_term_summary.iloc[:, reordered_indices]'''
            if z == 1:
                sum_1 = go_term_summary.sum(axis=1)
                go_term_summary1 = go_term_summary.div(go_term_summary.sum(axis=1), axis=0)
            elif z == 2:
                sum_2 = go_term_summary.sum(axis=1)
                go_term_summary2 = go_term_summary.div(go_term_summary.sum(axis=1), axis=0)

        # Combine the two dataframes by concatenating them side by side
        combined_go_term_summary = pd.concat([go_term_summary1, go_term_summary2], axis=1)

        sum_combined = pd.concat([sum_1, sum_2], axis=1)
        print(sum_combined)
        # append the sum of each row to the index
        new_index = []
        sum_combined.fillna(0, inplace=True)
        for index in combined_go_term_summary.index:
            if sum_combined.loc[index].iloc[0] > 0 and sum_combined.loc[index].iloc[1] > 0:
                new_index.append(f"{index} (n(cleavable) = {sum_combined.loc[index].iloc[0]}), n(noncleavable) = {sum_combined.loc[index].iloc[1]})")
            elif sum_combined.loc[index].iloc[0] > 0 and sum_combined.loc[index].iloc[1] == 0:
                new_index.append(f"{index} (n(cleavable) = {sum_combined.loc[index].iloc[0]})")
            elif sum_combined.loc[index].iloc[1] > 0 and sum_combined.loc[index].iloc[0] == 0:
                new_index.append(f"{index} (n(noncleavable) = {sum_combined.loc[index].iloc[1]})")
        combined_go_term_summary.index = new_index

        # plot the combined heatmap
        plt.figure(figsize=(18, 8))  # Increase width from 12 to 18 (or more if needed)
        norm = plt.Normalize(0, 1)
        ax = sns.heatmap(combined_go_term_summary, cmap="YlGnBu", annot=False, fmt=".2f", cbar=True, norm=norm)
        plt.title(f"GO Term Heatmap for {name} ({function})")
        #plt.xlabel("╚════════════════cleavable MTS════════════════╝╚═══════════════no cleavable MTS═══════════════╝")
        plt.xlabel(" cleavable MTS (left) vs non-cleavable MTS (right)")
        plt.ylabel("GO Term")
        plt.xticks(rotation=0, ha='right')  # Rotate and right-align x labels
        plt.yticks(rotation=0)
        plt.tight_layout()  # Adjust layout to prevent clipping
        plt.savefig(os.path.join(output_dir_per_organism, f"{name}_GO_term_heatmap_all_mito_split_into_cleavable_noncleavable.pdf"), dpi=300, bbox_inches='tight')
        plt.close()
        print(f"Heatmap saved for {name} with filter {filter_by}.")
        print(f"Processed {name} with filter {filter_by}.")
    
if __name__ == "__main__":
    # Define the names of the organisms
    organism_names = [
    "Homo_sapiens", "Mus_musculus",
    "Caenorhabditis_elegans", "Drosophila_Melanogaster", "Arabidopsis_thaliana""Saccharomyces_cerevisiae"]
    output_dir = "pipeline/output/output_20250513_133508"
    input_dir = "pipeline/input"
    filter = True
    absolute_threshold = 5
    relative_threshold = 0.5  # Minimum number of occurrences for a GO term to be included in the heatmap
    filter_by = "C" 
    run(organism_names, input_dir, output_dir, relative_threshold, filter, filter_by, absolute_threshold)
