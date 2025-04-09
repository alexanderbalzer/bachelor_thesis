import subprocess

def run_perl_script(perl_script_path, input_file, flag, output_file):
    """
    Run the MitoFates Perl script.
    """
    with open(output_file, "w") as output_handle:
        subprocess.run(
            ["perl", perl_script_path, input_file, flag],
            stdout=output_handle,
            stderr=subprocess.PIPE,
            check=True
        )

def parse_mts_annotations(mts_file):
    """
    Parse MTS-cleavable annotations from the MitoFates output.
    """
    mts_annotations = {}
    with open(mts_file, "r") as file:
        for line in file:
            parts = line.strip().split("\t")
            protein_id = parts[0]
            probability = float(parts[1])
            mts_annotations[protein_id] = probability
    return mts_annotations

def filter_proteins_by_mts(proteome, mts_annotations, threshold=0.9):
    """
    Filter proteins by MTS-cleavable probability.
    """
    filtered_proteins = []
    for protein_id, sequence in proteome:
        if protein_id in mts_annotations and mts_annotations[protein_id] >= threshold:
            filtered_proteins.append((protein_id, sequence))
    return filtered_proteins