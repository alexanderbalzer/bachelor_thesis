import os

def save_fasta(proteins, output_file):
    """
    Save proteins to a FASTA file.
    """
    with open(output_file, "w") as file:
        for protein_id, sequence in proteins:
            file.write(f">{protein_id}\n{sequence}\n")

def log_message(message):
    """
    Log a message to the console.
    """
    print(f"[LOG] {message}")