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

def transform_labels_to_names(organism_names):
    """
    Transform organism labels to names 
    """
    names = []
    for label in organism_names:
        # Replace underscores with spaces
        name = label.replace("_", " ")
        # Capitalize the first letter of each word
        name = name.title()
        # Add the transformed name to the list
        names.append(name)

    return names