from Bio import SeqIO
import pandas as pd
import pandas as pd
from Bio import SeqIO
import os
import logging
import subprocess





# This script filters protein sequences based on Gene Ontology (GO) terms and mitochondrial targeting sequences (MTS).
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

def parse_mts_cleavable_annotations(mts_file):
    """
    Parse MTS-cleavable annotations from a file.
    """
    cleavable = {}
    with open(mts_file, "r") as file:
      #  file_without_header = file.readlines()[2:]  
        for line in file:
            fields = line.strip().split("\t")
            if len(fields) < 3:
                continue
            protein_id = fields[0]
            cleavable[protein_id] = fields[2]
    logging.info(f"Parsed {len(cleavable)} MTS-cleavable annotations.")
    return cleavable

def parse_probability_of_mts (mts_file):
    """
    Parse probability of MTS-cleavable annotations from a file.
    """
    probability = {}
    with open(mts_file, "r") as file:
        file_without_header = file.readlines()[1:]  
        for line in file_without_header:
            fields = line.strip().split("\t")
            if len(fields) < 3:
                continue
            protein_id = fields[0]
            probability[protein_id] = float(fields[1])
    logging.info(f"Parsed {len(probability)} probability of MTS-cleavable annotations.")
    return probability


def add_mts_cleavable_to_dataframe(df, probability_of_mts):
    """
    Add MTS-cleavable information to the DataFrame based on the protein IDs.
    """
    for index, row in df.iterrows():
        protein_id = row["Header"]  
        if probability_of_mts[protein_id] > 0.9:
            df.at[index, "MTS-cleavable?"] = "Yes"
        elif probability_of_mts[protein_id] < 0.9:
            df.at[index, "MTS-cleavable?"] = "No"
        else:
            df.at[index, "MTS-cleavable?"] = "Unknown"
    return df

def filter_proteins_by_mts(dataframe, cleavable):
    cleavable_MTS = []
    for index, row in dataframe.iterrows():
        if row["MTS-cleavable?"] == cleavable:
            cleavable_MTS.append((row["Header"], row["Sequence"]))
    return cleavable_MTS

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

def run_perl_script(perl_script_path, input_file, flag, output_file):
    with open(output_file, "w") as output_file:
        subprocess.run(
            ["perl", perl_script_path, input_file, flag],
            stdout=output_file,
            stderr=subprocess.PIPE,
            check=True
    )


logging.basicConfig(level=logging.INFO)


logging.info("Starting the pipeline...")

# Define the organisms to process
organisms = ["human"]
cleavable = "Yes" #Yes or No
threshold = 0.9 #threshold for probability of MTS-cleavable
flag = "metazoa"

perl_script_path = "/home/abalzer/Downloads/MitoFates_1.2/MitoFates/MitoFates.pl"



for i, name in enumerate(organisms, start=1):
    # Log the current organism being processed
    logging.info(f"Processing organism: {name}")
    # Log the current iteration
    logging.info(f"Durchlauf[{i}/{len(organisms)}]")
    # Define the input files
    
    annotation_file = "pipeline/input/" + str(name) + ".goa"  
    go_annotation = parse_go_annotations(annotation_file)

    fasta_file = "pipeline/input/" + str(name) + ".fasta"
        
    # Define the output file
    output_file = f"output files/filtered_proteins_by_GO_and_mts_subset_for_{name}.fasta"


    proteome = fasta_to_dataframe(fasta_file)

    proteome_with_go_terms = add_go_terms_to_dataframe(proteome, go_annotation)
    #print(proteome_with_go_terms)

    target_go_term = "GO:0005739" #go term for mitochondrion

    filtered_proteins = filter_proteins_by_go(proteome_with_go_terms, target_go_term)
    #print(filtered_proteins)
    #filtered_proteins = filtered_proteins[:2000]  # Limit to the first 2000 proteins

    valid_amino_acids = set("ACDEFGHIKLMNPQRSTVWY") #ARNDCEQGHILKMFPSTWYV


    # Filter out proteins with invalid amino acids
    filtered_proteins = [
        (protein_id, protein_seq)
        for protein_id, protein_seq in filtered_proteins
        if set(protein_seq).issubset(valid_amino_acids)
    ]


    with open("pipeline/cache/filtered_proteins_by_GO_for_" + str(name) + ".fasta", "w") as output_handle:
        for protein_id, protein_seq in filtered_proteins:
            output_handle.write(f">{protein_id}\n{protein_seq}\n")
    logging.info(f"Filtered proteins by GO terms and saved to {output_file}")

    input_file = "pipeline/cache/filtered_proteins_by_GO_for_" + str(name) + ".fasta"
    output_file = "pipeline/cache/mito_fates_for_" + str(name) + ".cgi"

    logging.info(f"Running MitoFates Perl script for {name}...")
    run_perl_script(perl_script_path, input_file, flag, output_file)


    fasta_file = "pipeline/cache/filtered_proteins_by_GO_for_" + str(name) + ".fasta"
    proteome = fasta_to_dataframe(fasta_file)

    mts_cleavable_file = "pipeline/cache/mito_fates_for_" + str(name) + ".cgi"

    check_file_exists(fasta_file)
    check_file_exists(mts_cleavable_file)

    mts_cleavable = parse_mts_cleavable_annotations(mts_cleavable_file)
    logging.info(f"Parsed {len(mts_cleavable)} MTS-cleavable annotations.")
    probability_of_mts = parse_probability_of_mts(mts_cleavable_file)
    logging.info(f"Parsed {len(probability_of_mts)} probability of MTS-cleavable annotations.")
    proteome = add_mts_cleavable_to_dataframe(proteome, probability_of_mts)

    filtered_proteins = filter_proteins_by_mts(proteome, cleavable)


    if cleavable == "No":
        output_file = "pipeline/output/filtered_by_GO_no_cleavable_mts_for_" + str(name) + ".fasta"
    if cleavable == "Yes":
        output_file = "pipeline/output/filtered_by_GO_cleavable_mts_for_" + str(name) + ".fasta"

    with open(output_file, "w") as f:
        for header, sequence in filtered_proteins:
            f.write(f">{header}\n{sequence}\n")
    logging.info(f"Filtered proteins by MTS-cleavable status and saved to {output_file}")

    os.remove("pipeline/cache/filtered_proteins_by_GO_for_" + str(name) + ".fasta")
    os.remove("pipeline/cache/mito_fates_for_" + str(name) + ".cgi")
    logging.info(f"Removed temporary files for {name}")
logging.info("Pipeline completed successfully.")










