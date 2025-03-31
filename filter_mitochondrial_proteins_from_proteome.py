from Bio import SeqIO
import pandas as pd


def parse_go_annotations(annotation_file):
    """
    Parse a GO annotation file (e.g., GAF format) to create a dictionary
    mapping protein IDs to their associated GO terms.
    """
    go_annotation = {}
    with open(annotation_file, "r") as file:
        for line in file:
            if line.startswith("!"):  
                continue
            fields = line.strip().split("\t")
            protein_id = fields[1]  #protein ID (column 2 in GAF)
            function = fields[8]  #function type (column 9 in GAF)
            if function != "C":  #filter for location information
                continue
            go_term = fields[4]     #GO term (column 5 in GAF)
            if protein_id not in go_annotation:
                go_annotation[protein_id] = []
            go_annotation[protein_id].append(go_term)
    return go_annotation


def fasta_to_dataframe(fasta_file):
    headers = []
    annotations = []
    sequences = []

    for record in SeqIO.parse(fasta_file, "fasta"):
        headers.append(record.id)
        sequences.append(str(record.seq))
    
    df = pd.DataFrame({
        "Header": headers,
        "GO_Term": [""] * len(headers),    
        "Sequence": sequences
    })
    
    return df

def add_go_terms_to_dataframe(df, go_annotation):
    """
    Add GO terms to the DataFrame based on the protein IDs.
    """
    for index, row in df.iterrows():
        protein_id = row["Header"].strip().split("|")[1]  
        if protein_id in go_annotation:
            go_terms = go_annotation[protein_id]
            df.at[index, "GO_Term"] = ", ".join(go_terms)
        else:
            df.at[index, "GO_Term"] = "No GO term found"
    return df

def filter_proteins_by_go(dataframe, target_go_term):
    filtered_proteins = []
    for index, row in dataframe.iterrows():
        go_terms = row["GO_Term"]
        if target_go_term in go_terms:
            filtered_proteins.append((row["Header"], row["Sequence"]))

    return filtered_proteins

annotation_file = "input files/9.C_elegans.goa"  
go_annotation = parse_go_annotations(annotation_file)

fasta_file = "input files/uniprotkb_proteome_UP000001940_AND_prot_2025_03_28.fasta"

proteome = fasta_to_dataframe(fasta_file)

proteome_with_go_terms = add_go_terms_to_dataframe(proteome, go_annotation)
print(proteome_with_go_terms)

target_go_term = "GO:0005739"

filtered_proteins = filter_proteins_by_go(proteome_with_go_terms, target_go_term)
print(filtered_proteins)

'''
with open("filtered_proteins.fasta", "w") as output_handle:
    for protein_id, protein_seq in filtered_proteins:
        output_handle.write(f">{protein_id}\n{protein_seq}\n")'
'''
