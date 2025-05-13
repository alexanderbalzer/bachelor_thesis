import requests
from Bio import SeqIO
import pandas as pd
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

proteins = [ "RPLP0", "RPL3", "RPL4", "RPL5", "RPL6", "RPL7A", "RPL7", "RPL8", "RPL9", "RPL10", 
"RPL11", "RPL13A", "RPL13", "RPL14", "RPL15", "RPL17", 
"RPL18A", "RPL18", "RPL19", "RPL22", "RPL23A", "RPL23", "RPL24", "RPL26", "RPL27A",
"RPL27", "RPL28", "RPL29", "RPL30", "RPL31", "RPL32", "RPL34", "RPL35A", "RPL35", "RPL36A",
"RPL36", "RPL37A", "RPL37", "RPL38", "RPL39", "RPL40" ]

def fetch_protein_sequences(proteins, fasta_file):
    protein_sequences = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        protein_description = record.description
        informations = protein_description.split(" ")
        for info in informations:
            if info.startswith("GN="):
                gene_id = info.split("=")[1]
                if gene_id in proteins:
                    protein_sequences[gene_id] = str(record.seq)
                    break
    return protein_sequences

fasta_file_path = "pipeline/input/Homo_sapiens.fasta"  # Replace with the actual path to your proteome fasta file
protein_sequences = fetch_protein_sequences(proteins, fasta_file_path)
# convert to panda
df = pd.DataFrame(list(protein_sequences.items()), columns=["gene_id", "sequence"])

# Convert DataFrame to list of SeqRecord objects
records = [
    SeqRecord(Seq(row["sequence"]), id=row["gene_id"], description="")
    for _, row in df.iterrows()
]

# Save as FASTA
with open("protein_sequences_60s.fasta", "w") as fasta_out:
    SeqIO.write(records, fasta_out, "fasta")