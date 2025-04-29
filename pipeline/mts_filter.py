import subprocess
import os
import pandas as pd
from Bio import SeqIO
import logging
import shutil

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





def run(list_of_organisms, output_dir, cleavable, mitofates_path, flagdict, delete_cache, threshold,run_from_scratch, amount_of_proteins_per_step, last_run):
    """
    Run the MitoFates Perl script and filter proteins by MTS-cleavable probability.
    """
    # Check if the Perl script exists
    check_file_exists(mitofates_path)
    logging.info(f"Perl script found at {mitofates_path}")

    # Run the Perl script for each organism
    for organism in list_of_organisms:
        output_dir_per_organism = output_dir + "/" + organism 
        logging.info(f"Running MitoFates for {organism}")
        if not run_from_scratch:    
            if os.path.exists(os.path.join(last_run, f"{organism}_filtered_by_go_and_mts.fasta")):
                # copy the file to the new cache directory
                shutil.copy(os.path.join(last_run, f"{organism}_filtered_by_go_and_mts.fasta"), os.path.join(output_dir_per_organism, f"{organism}_filtered_by_go_and_mts.fasta"))
                logging.info(f"File already exists for {organism}: {os.path.join(last_run, f'{organism}_filtered_by_go_and_mts.fasta')}")
                continue

            # check if the output file already exists and skip if it does
            output_file = os.path.join(output_dir_per_organism, f"{organism}_filtered_by_go_and_mts.fasta")
            if os.path.exists(output_file):
                logging.info(f"Output file already exists for {organism}: {output_file}")
                continue
        try:
            input_file = os.path.join(output_dir_per_organism, f"filtered_proteins_by_GO_for_{organism}.fasta")
            output_mitofates__file = os.path.join(output_dir_per_organism, f"mitofates_for_{organism}.cgi")
            run_perl_script(mitofates_path, input_file, flagdict[organism], output_mitofates__file)
            logging.info(f"Perl script completed for {organism}")
        except subprocess.CalledProcessError as e:
            logging.error(f"Error running MitoFates Perl script for {organism}: {e.stderr.decode()}")
            continue
        except FileNotFoundError:
            logging.error(f"Input file not found for {organism}: {input_file}")
            continue
        proteome = fasta_to_dataframe(input_file)
        probability_of_mts = parse_probability_of_mts(output_mitofates__file)
        proteome = add_mts_cleavable_to_dataframe(proteome, probability_of_mts, threshold)
        filtered_proteins_cleavable = filter_proteins_by_mts(proteome, cleavable="Yes")
        filtered_proteins_non_cleavable = filter_proteins_by_mts(proteome, cleavable="No")
        amount_of_proteins_per_step.at["Mitochondrial with MTS", organism] = len(filtered_proteins_cleavable)
        with open(os.path.join(output_dir_per_organism, f"/{organism}_filtered_by_GO_cleavable_mts.fasta"), "w") as output_handle:
            for header, sequence in filtered_proteins_cleavable:
                output_handle.write(f">{header}\n{sequence}\n")
        logging.info(f"Filtered proteins saved to {output_dir_per_organism}/{organism}_filtered_by_GO_cleavable_mts.fasta")
        with open(os.path.join(output_dir_per_organism, f"/filtered_proteins_by_GO_noncleavable_mts_{organism}.fasta"), "w") as output_handle:
            for header, sequence in filtered_proteins_non_cleavable:
                output_handle.write(f">{header}\n{sequence}\n")

        if delete_cache == "yes":
            os.remove(input_file)
            os.remove(output_file)
            os.remove(os.path.join(output_dir_per_organism, f"filtered_proteins_by_GO_for_{organism}.fasta"))
            logging.info(f"Deleted cache files for {organism}")

    logging.info("MitoFates filtering completed.")
    return amount_of_proteins_per_step

if __name__ == "__main__":
    # Example usage
    organism_names = ["Arabidopsis_thaliana",
    "Caenorhabditis_elegans",
    "Candida_glabrata",
    "Clavispora_lusitaniae",
    "Debaryomyces_hansenii",
    "Drosophila_Melanogaster",
    "Geotrichum_candidum",
    "human",
    "human_with_isoforms",
    "Lachancea_thermotolerans",
    "Mus_musculus",
    "Physcomitrium_patens",
    "Saccharomyces_cerevisiae",
    "Scheffersomyces_stipitis",
    "Schizosaccharomyces_pombe",
    "Yarrowia_lipolytica",
    "Zygosaccharomyces_rouxii"]
    input_dir = "pipeline/input"
    output_dir_per_organism = "pipeline/cache/cache_20250414_224234/"
    output_dir = "pipeline/output/output_20250414_224234"
    create_heatmap = True
    heatmap_type = "hgt"
    create_phylogenetic_tree = True
    phylo_tree_type = "hgt"
    reference = "proteome"
    cleavable = "Yes"
    perl_script_path = "/home/abalzer/MitoFates/MitoFates.pl"
    with open("pipeline/flaglist.txt", "r") as file:
        flaglist = {}
        for line in file:
            key, value = line.strip().split(":")
            flaglist[key.strip()] = value.strip()
    delete_cache = "No"
    threshold = 0.9
    run_from_scratch = True
    amount_of_proteins_per_step = pd.DataFrame(index=["Start", "Mitochondrial", "Mitochondrial with MTS"], columns=organism_names)
    last_run = "pipeline/cache/cache_20250414_224234/"
    run(organism_names, output_dir, cleavable, perl_script_path, flaglist, delete_cache, threshold, run_from_scratch, amount_of_proteins_per_step, last_run)