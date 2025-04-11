import subprocess
import os
import pandas as pd
from Bio import SeqIO
import logging

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
    return probability


def add_mts_cleavable_to_dataframe(df, probability_of_mts, threshold):
    for index, row in df.iterrows():
        protein_id = row["Header"]
        probability = probability_of_mts.get(protein_id, 0.0)  # Default to 0.0 if not found
        if probability >= threshold:
            df.at[index, "MTS-cleavable?"] = "Yes"
        elif probability < threshold:
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
        logging.error(f"File not found: {file_path}")
        raise

def run_perl_script(mitofates_path, input_file, flag, output_file):

    env = os.environ.copy()
    env["PATH"] = "/usr/bin:" + env["PATH"]  
    with open(output_file, "w") as output_handle:
        subprocess.run(
            ["perl", mitofates_path, input_file, flag],
            stdout=output_handle,
            stderr=subprocess.PIPE,
            env=env,
            check=True
        )
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





def run(list_of_organisms, cache_dir, output_dir, cleavable, mitofates_path, flagdict, delete_cache, threshold):
    """
    Run the MitoFates Perl script and filter proteins by MTS-cleavable probability.
    """
    # Check if the Perl script exists
    check_file_exists(mitofates_path)
    logging.info(f"Perl script found at {mitofates_path}")

    # Run the Perl script for each organism
    for organism in list_of_organisms:
        # check if the output file already exists and skip if it does
        output_file = os.path.join(output_dir, f"{organism}_filtered_by_go_and_mts.fasta")
        if os.path.exists(output_file):
            logging.info(f"Output file already exists for {organism}: {output_file}")
            continue
        try:
            input_file = os.path.join(cache_dir, f"filtered_proteins_by_GO_for_{organism}.fasta")
            output_file = os.path.join(cache_dir, f"mitofates_for_{organism}.cgi")
            run_perl_script(mitofates_path, input_file, flagdict[organism], output_file)
            logging.info(f"Perl script completed for {organism}")
        except subprocess.CalledProcessError as e:
            logging.error(f"Error running MitoFates Perl script for {organism}: {e.stderr.decode()}")
            continue
        except FileNotFoundError:
            logging.error(f"Input file not found for {organism}: {input_file}")
            continue
        proteome = fasta_to_dataframe(input_file)
        probability_of_mts = parse_probability_of_mts(output_file)
        proteome = add_mts_cleavable_to_dataframe(proteome, probability_of_mts, threshold)
        filtered_proteins = filter_proteins_by_mts(proteome, cleavable)

        with open(os.path.join(output_dir, f"{organism}_filtered_by_go_and_mts.fasta"), "w") as output_handle:
            for header, sequence in filtered_proteins:
                output_handle.write(f">{header}\n{sequence}\n")
        logging.info(f"Filtered proteins saved to {output_dir}{organism}_filtered_by_go_and_mts.fasta")

        if delete_cache == "yes":
            os.remove(input_file)
            os.remove(output_file)
            os.remove(os.path.join(cache_dir, f"filtered_proteins_by_GO_for_{organism}.fasta"))
            logging.info(f"Deleted cache files for {organism}")

    logging.info("MitoFates filtering completed.")

