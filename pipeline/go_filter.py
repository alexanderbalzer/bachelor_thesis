from Bio import SeqIO
import pandas as pd
import pandas as pd
from Bio import SeqIO
import os
import logging
import shutil


def parse_go_annotations(annotation_file):
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
#            if function != "C":  #filter for location information
#                continue
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
        protein_id = row["Header"].strip().split("|")[1] 
        if protein_id in go_annotation:
            go_terms = go_annotation[protein_id]
            df.at[index, "GO_Term"] = ", ".join(go_terms)
        else:
            df.at[index, "GO_Term"] = "No GO term found"
    return df

def filter_proteins_by_go(dataframe, target_go_term, not_in_go_term):
    logging.info(f"Filtering proteins by GO term: {target_go_term}")
    filtered_proteins = []
    print(f"Filtering proteins by GO term: {target_go_term}")
    print(f"leave out GO term: {not_in_go_term}")
    for index, row in dataframe.iterrows():
        if not_in_go_term is not False:
            go_terms = str(row["GO_Term"]).split(", ")
            if str(target_go_term) in go_terms and str(not_in_go_term) not in go_terms:
                filtered_proteins.append((row["Header"], row["Sequence"]))
        else:
            go_terms = str(row["GO_Term"]).split(", ")
            if str(target_go_term) in go_terms:
                filtered_proteins.append((row["Header"], row["Sequence"]))
    return filtered_proteins

def fasta_to_dataframe(fasta_file):
    headers = []
    sequences = []

    for record in SeqIO.parse(fasta_file, "fasta"):
        headers.append(record.id)
        sequences.append(str(record.seq))
    
    df = pd.DataFrame({
        "Header": headers,
        "threshold": [""] * len(headers),
        "MTS-cleavable?": [""] * len(headers),
        "Sequence": sequences
    })
    return df
def check_file_exists(file_path):
    """
    Check if a file exists.
    """
    try:
        with open(file_path, 'r') as file:
            pass
    except FileNotFoundError:
        print(f"File not found: {file_path}")
        raise


def run(list_of_organisms, input_dir, output_dir, target_go_term, run_from_scratch, amount_of_proteins_per_step, last_run, not_in_go_term):
    """
    Main function to run the pipeline.
    """
    for i, name in enumerate(list_of_organisms, start=1):
        logging.info(f"Filtering by GO term for: {name}")
        logging.info(f"Run[{i}/{len(list_of_organisms)}]")

        # Define the input and output files
        annotation_file = os.path.join(input_dir, f"{name}.goa")
        fasta_file = os.path.join(input_dir, f"{name}.fasta")
        output_filtered_by_GO_file = output_dir + "/" + name + f"/filtered_proteins_by_GO_for_{name}.fasta"
        # Create the output directory if it doesn't exist
        os.makedirs(os.path.dirname(output_filtered_by_GO_file), exist_ok=True)

        last_run_output_file = last_run + "/" + name + f"/filtered_proteins_by_GO_for_{name}.fasta"

        if not run_from_scratch:
        # check if the output file already exists
            if os.path.exists(last_run_output_file):
                # copy the file to the new cache directory
                shutil.copy(last_run_output_file, output_filtered_by_GO_file)
                logging.info(f"Output file {output_filtered_by_GO_file} already exists. Skipping.")
                continue

        # Check if the input files exist
        check_file_exists(annotation_file)

        # Parse the GO annotations and the FASTA file
        logging.info(f"Parsing GO annotations from {annotation_file}...")
        go_annotation = parse_go_annotations(annotation_file)
        if not go_annotation:
            logging.error(f"No GO annotations found for {name}. Check the .goa file.")
            continue
        proteome = fasta_to_dataframe(fasta_file)
        amount_of_proteins_per_step.at["Start", name] = len(proteome)

        proteome_with_go_terms = add_go_terms_to_dataframe(proteome, go_annotation)
        if proteome_with_go_terms.empty:
            logging.error(f"No proteins found in the FASTA file {fasta_file}.")
            continue

        filtered_proteins = filter_proteins_by_go(proteome_with_go_terms, target_go_term, not_in_go_term)
        if not filtered_proteins:
            logging.error(f"No proteins matched the target GO term {target_go_term} for {name}.")

        valid_amino_acids = set("ACDEFGHIKLMNPQRSTVWY") #ARNDCEQGHILKMFPSTWYV


        # Filter out proteins with invalid amino acids
        filtered_proteins = [
            (protein_id, protein_seq)
            for protein_id, protein_seq in filtered_proteins
            if set(protein_seq).issubset(valid_amino_acids)
        ]
        if not filtered_proteins:
            logging.error(f"All proteins for {name} were filtered out due to invalid amino acids.")

        amount_of_proteins_per_step.at["Mitochondrial", name] = len(filtered_proteins)
        # Write the filtered proteins to the cache
        with open(output_filtered_by_GO_file, "w") as output_handle:
            for protein_id, protein_seq in filtered_proteins:
                output_handle.write(f">{protein_id}\n{protein_seq}\n")
        logging.info(f"Filtered proteins by GO term and saved to {output_filtered_by_GO_file}")
    return amount_of_proteins_per_step






