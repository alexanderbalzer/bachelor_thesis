from Bio import SeqIO

def parse_go_annotations(annotation_file):
    """
    Parse GO annotations from a file.
    """
    go_annotations = {}
    with open(annotation_file, "r") as file:
        for line in file:
            parts = line.strip().split("\t")
            protein_id = parts[0]
            go_terms = parts[1:]
            go_annotations[protein_id] = go_terms
    return go_annotations

def filter_proteins_by_go(proteome, go_annotations, target_go_term):
    """
    Filter proteins by the target GO term.
    """
    filtered_proteins = []
    for protein_id, sequence in proteome:
        if protein_id in go_annotations and target_go_term in go_annotations[protein_id]:
            filtered_proteins.append((protein_id, sequence))
    return filtered_proteins

def load_proteome(fasta_file):
    """
    Load proteome sequences from a FASTA file.
    """
    proteome = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        proteome.append((record.id, str(record.seq)))
    return proteome