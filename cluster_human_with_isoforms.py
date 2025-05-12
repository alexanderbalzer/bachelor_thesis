from Bio import SeqIO
import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np

def load_fasta_to_dataframe(fasta_file):
    # Parse the FASTA file and extract protein ID and first 50 amino acids
    data = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        protein_id = record.id.strip().split("|")[1]  # Extract the protein ID
        first_aa = str(record.seq[:length])
        data.append((protein_id, first_aa))
    
    # Create a pandas DataFrame
    df = pd.DataFrame(data, columns=["Protein_ID", "First_AA"])
    return df



if __name__ == "__main__":
    length = 70  # Length of the sequence to consider
    # Path to the FASTA file
    fasta_file = "pipeline/input/human_with_isoforms.fasta"
    # Load the FASTA file into a DataFrame
    df = load_fasta_to_dataframe(fasta_file)
    before_group = len(df)

    # group by the first 50 amino acids and count occurrences, but add the protein IDs
    grouped_df = df.groupby("First_AA")["Protein_ID"].apply(list).reset_index()
    grouped_df["Count"] = grouped_df["Protein_ID"].apply(len)
    grouped_df = grouped_df.sort_values(by="Count", ascending=False)
    print(grouped_df.head(10))
    
    after_group = len(grouped_df)
    print(f"Number of unique sequences before grouping: {before_group}")
    print(f"Number of unique sequences after grouping: {after_group}")
    # group by the count of sequences
    occurences = grouped_df.groupby("Count").size().reset_index(name="Frequency")
    # Plot the histogram of occurrences
    plt.figure(figsize=(10, 6))
    plt.bar(occurences["Count"], occurences["Frequency"], color='blue', alpha=0.7, edgecolor='black')
    plt.yscale('log')
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.tight_layout()
    plt.xlabel("Isoforms per Gene")
    plt.ylabel("Frequency")
    plt.title("Histogram of Sequence Counts")
    
    plt.show()

    # Save the grouped DataFrame to a CSV file
    '''output_file = "pipeline/input/human_with_isoforms_grouped.csv"
    grouped_df.to_csv(output_file, index=False)
    print(f"Grouped sequences saved to {output_file}")'''
