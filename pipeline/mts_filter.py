import subprocess
import os
import pandas as pd

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


def add_mts_cleavable_to_dataframe(df, probability_of_mts):
    for index, row in df.iterrows():
        protein_id = row["Header"]
        probability = probability_of_mts.get(protein_id, 0.0)  # Default to 0.0 if not found
        if probability > 0.9:
            df.at[index, "MTS-cleavable?"] = "Yes"
        elif probability < 0.9:
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
    with open(output_file, "w") as output_handle:
        subprocess.run(
            ["perl", perl_script_path, input_file, flag],
            stdout=output_handle,
            stderr=subprocess.PIPE,
            check=True
        )






def run(list_of_organisms, input_dir, cache_dir, output_dir, cleavable, mitofates_path, flaglist):
    """
    Run the MitoFates Perl script and filter proteins by MTS-cleavable probability.
    """
    # Check if the Perl script exists
    check_file_exists(mitofates_path)

    # Run the Perl script for each organism
    for organism in list_of_organisms:
        input_file = os.path.join(input_dir, f"{organism}.fasta")
        output_file = os.path.join(cache_dir, f"mitofates_{organism}.txt")
        run_perl_script(mitofates_path, input_file, flaglist[organism], output_file)
        print(f"Perl script completed for {organism}")

    # Parse the MTS-cleavable annotations
    cleavable_annotations = {}
    for organism in list_of_organisms:
        mts_file = os.path.join(cache_dir, f"mitofates_{organism}.txt")
        cleavable_annotations[organism] = parse_mts_cleavable_annotations(mts_file)












