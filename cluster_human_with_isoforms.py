from Bio import SeqIO
import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np
import tqdm as tqdm
from tqdm import tqdm

def load_fasta_to_dataframe(fasta_file, length):
    # Parse the FASTA file and extract protein ID and first 50 amino acids
    data = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        protein_id = record.id  # Extract the protein ID
        first_aa = str(record.seq)
        sequence = str(record.seq)
        protein_evidence = None
        if "PE=" in record.description:
            try:
                protein_evidence = int(record.description.split("PE=")[1].split()[0])
            except ValueError:
                protein_evidence = None
        data.append((protein_id, first_aa, sequence, protein_evidence))
    
    # Create a pandas DataFrame
    df = pd.DataFrame(data, columns=["Protein_ID", "First_AAs", "sequence", "protein_evidence"])
    return df


def run(organism_names, output_dir):
    comparison_before_after = []
    amount_of_proteins_organisms = []
    for name in organism_names:
        length = 60  # Length of the sequence to consider
        # Path to the FASTA file
        fasta_file = f"pipeline/input/{name}.fasta"
        # Load the FASTA file into a DataFrame
        df = load_fasta_to_dataframe(fasta_file, length)
        # set the index to the protein ID
        df.set_index("Protein_ID", inplace=True)
        df["Protein_ID"] = df.index
        before_group = len(df)

        grouped_df = df.copy()
        after_group = len(grouped_df)
        # Filter out rows with protein_evidence smaller than 4
        grouped_df5 = grouped_df[grouped_df["protein_evidence"] <= 5]
        grouped_df4 = grouped_df[grouped_df["protein_evidence"] <= 4]
        grouped_df3 = grouped_df[grouped_df["protein_evidence"] <= 3]
        grouped_df2 = grouped_df[grouped_df["protein_evidence"] <= 2]
        grouped_df1 = grouped_df[grouped_df["protein_evidence"] <= 1]
        
        
        # group by the first 60 amino acids and count occurrences, but add the protein IDs
        grouped_df = grouped_df1.groupby("First_AAs")["Protein_ID"].apply(list).reset_index()
        grouped_df["Count"] = grouped_df["Protein_ID"].apply(len)
        #grouped_df = grouped_df.sort_values(by="Count", ascending=False)
        # group by the count of sequences
        occurences = grouped_df.groupby("Count").size().reset_index(name="Frequency")
        # Plot the histogram of occurrences
        plt.figure(figsize=(10, 6))
        plt.bar(occurences["Count"], occurences["Frequency"], color="blue", alpha=0.7)
        plt.yscale('log')
        plt.grid(axis='y', linestyle='--', alpha=0.7)
        plt.xlabel("Isoforms per Gene")
        plt.ylabel("Frequency")

        # Add the sequence of the first Protein_ID in each group to the grouped_df
        grouped_df["sequence"] = grouped_df["Protein_ID"].apply(lambda ids: df.loc[ids[0], "sequence"])
        grouped_df["Protein_ID"] = grouped_df["Protein_ID"].apply(lambda ids: ids[0])

        amount_of_proteins = {
            "organism": name,
            "5": len(grouped_df5),
            "4": len(grouped_df4),
            "3": len(grouped_df3),
            "2": len(grouped_df2),
            "1": len(grouped_df1),
            "unique_sequences": len(grouped_df),
            "amount_of_duplicates": len(grouped_df1) - len(grouped_df)}
        amount_of_proteins_organisms.append(amount_of_proteins)

        os.makedirs(output_dir, exist_ok=True)

        # Save the original sequences to a fasta file
        original_fasta_folder = os.path.join(output_dir, "original_sequences")
        os.makedirs(original_fasta_folder, exist_ok=True)
        original_fasta_file = os.path.join(original_fasta_folder, f"{name}_original.fasta")
        with open(original_fasta_file, "w") as f:
            for index, row in df.iterrows():
                f.write(f">{row['Protein_ID']}\n")
                f.write(f"{row['sequence']}\n")
        # Save the grouped DataFrame to a fasta file
        grouped_fasta_folder = os.path.join(output_dir, "grouped_sequences")
        os.makedirs(grouped_fasta_folder, exist_ok=True)
        grouped_fasta_file = os.path.join(grouped_fasta_folder, f"{name}.fasta")
        with open(grouped_fasta_file, "w") as f:
            for index, row in grouped_df.iterrows():
                f.write(f">{row['Protein_ID']}\n")
                f.write(f"{row['sequence']}\n")
    

        # save the histogram to a file
        histogram_folder = os.path.join(output_dir, "histograms")
        os.makedirs(histogram_folder, exist_ok=True)
        # Save the histogram as a PDF file
        histogram_file = os.path.join(histogram_folder, f"{name}_histogram.pdf")
        plt.title(f"Histogram of Isoforms per Gene for {name}")
        plt.tight_layout()
        plt.savefig(histogram_file, dpi=300)
        plt.close() 
    amount_of_proteins_df = pd.DataFrame(amount_of_proteins_organisms)
    # Save the amount of proteins DataFrame to a CSV file
    amount_of_proteins_file = os.path.join(output_dir, "amount_of_proteins.csv")
    amount_of_proteins_df.to_csv(amount_of_proteins_file, index=False)

if __name__ == "__main__":
    # Example usage
    output_dir = "pipeline/grouping_of_proteins_with_same_n_terminus"
    organism_names = [
    "Homo_sapiens","Mus_musculus", "Rattus_norvegicus", "Dario_rerio",
    "Caenorhabditis_elegans", "Drosophila_Melanogaster", "Arabidopsis_thaliana", 
    "Saccharomyces_cerevisiae"]
    run(organism_names=organism_names, output_dir=output_dir)



