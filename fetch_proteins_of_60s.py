import requests
from Bio import SeqIO
import pandas as pd
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from io import StringIO

proteins = ["TOM5", "TOM6", "TOM7", "TOM20", "TOM22", "TOM40", "NAA10",
"NAA15", "HYPK", "HTT", "HSPA9", "VDAC2", "NACA", "BTF3"]

def fetch_protein_sequence(protein_name):
    url = f"https://rest.uniprot.org/uniprotkb/search?query={protein_name}&format=fasta&size=1"
    response = requests.get(url)
    if response.status_code == 200:
        return response.text
    else:
        print(f"Failed to fetch {protein_name}: {response.status_code}")
        return None

proteins = [protein.upper() + "_HUMAN" for protein in proteins]  # Ensure all protein names are uppercase
protein_sequences = []
for protein in proteins:
    sequence = fetch_protein_sequence(protein)
    if sequence:
        # Parse the FASTA format and create a SeqRecord
        handle = StringIO(sequence)  # ← FIX
        seq_record = SeqIO.read(handle, "fasta")
        protein_sequences.append(str(seq_record.seq))
# Save fetched sequences to a file
with open("fetched_proteins.fasta", "w") as fasta_file:
    for protein, seq in zip(proteins, protein_sequences):
        fasta_file.write(f"{seq}\n")