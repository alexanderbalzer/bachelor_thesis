from Bio import SeqIO
import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np

def load_fasta_to_dataframe(fasta_file):
    # Parse the FASTA file and extract protein ID and first 50 amino acids
    data = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        protein_id = record.id  # Extract the protein ID
        first_aa = str(record.seq[:length])
        sequence = str(record.seq)
        data.append((protein_id, first_aa, sequence))
    
    # Create a pandas DataFrame
    df = pd.DataFrame(data, columns=["Protein_ID", "First_AA", "sequence"])
    return df



if __name__ == "__main__":
    length = 60  # Length of the sequence to consider
    # Path to the FASTA file
    fasta_file = "pipeline/input/human_with_isoforms.fasta"
    # Load the FASTA file into a DataFrame
    df = load_fasta_to_dataframe(fasta_file)
    # set the index to the protein ID
    df.set_index("Protein_ID", inplace=True)
    df["Protein_ID"] = df.index
    before_group = len(df)

    # group by the first 50 amino acids and count occurrences, but add the protein IDs
    grouped_df = df.groupby("First_AA")["Protein_ID"].apply(list).reset_index()
    grouped_df["Count"] = grouped_df["Protein_ID"].apply(len)
    #grouped_df = grouped_df.sort_values(by="Count", ascending=False)
    print(grouped_df.head(10))
    
    after_group = len(grouped_df)
    print(f"Number of unique sequences before grouping: {before_group}")
    print(f"Number of unique sequences after grouping: {after_group}")
    # group by the count of sequences
    occurences = grouped_df.groupby("Count").size().reset_index(name="Frequency")
    # Plot the histogram of occurrences
    plt.figure(figsize=(10, 6))
    plt.hist(occurences["Count"], bins=5, weights=occurences["Frequency"], color="blue", alpha=0.7)
    plt.yscale('log')
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.xlabel("Isoforms per Gene")
    plt.ylabel("Frequency")

    # Add the sequence of the first Protein_ID in each group to the grouped_df
    grouped_df["sequence"] = grouped_df["Protein_ID"].apply(lambda ids: df.loc[ids[0], "sequence"])
    grouped_df["Protein_ID"] = grouped_df["Protein_ID"].apply(lambda ids: ids[0])
    print(grouped_df.head(10))

    # Save the grouped DataFrame to a fasta file
    output_dir = "pipeline/human_with_isoforms_grouping_infos"
    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(output_dir, "grouped_human_with_isoforms.fasta")
    with open(output_file, "w") as f:
        for index, row in grouped_df.iterrows():
            f.write(f">{row['Protein_ID']}\n")
            f.write(f"{row['sequence']}\n")
    
    # save the histogram to a file
    histogram_file = os.path.join(output_dir, "histogram.pdf")
    plt.savefig(histogram_file, dpi=300)
    plt.close() 


