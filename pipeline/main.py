import init
import os
from go_filter import parse_go_annotations, filter_proteins_by_go, load_proteome
from mts_filter import run_perl_script, parse_mts_annotations, filter_proteins_by_mts
from utils import save_fasta, log_message

def main():
    """
    Main function to run the pipeline.
    """
    # Check for config file and load it
    clearance = init.check_config()
    if not clearance:
        print("Configuration check failed. Exiting.")
        return
    
    # Load the config file
    config = init.load_config()
    input_dir = config['DEFAULT']['input_dir']
    output_dir = config['DEFAULT']['output_dir']
    cache_dir = config['DEFAULT']['cache_dir']
    target_go_term = config['DEFAULT']['target_GO_term']
    perl_script_path = config['DEFAULT']['mitofates_path']
    cleavable = config['DEFAULT']['cleavable']
    threshold = config['DEFAULT'].getfloat('threshold', 0.9)


    # Log the start of the pipeline
    log_message("Starting the pipeline...")

    # Read organism names from FASTA file names
    fasta_files = [f for f in os.listdir(input_dir) if f.endswith(".fasta")]
    organism_names = [os.path.splitext(f)[0] for f in fasta_files]
    log_message(f"Organism names extracted: {', '.join(organism_names)}")

    for organism in organism_names:
        log_message(f"Processing organism: {organism}")
        organism_input_dir = os.path.join(input_dir, organism)
        organism_output_dir = os.path.join(output_dir, organism)
        organism_cache_dir = os.path.join(cache_dir, organism)

        # Ensure organism-specific directories exist
        os.makedirs(organism_output_dir, exist_ok=True)
        os.makedirs(organism_cache_dir, exist_ok=True)

        # Update paths for organism-specific processing
        annotation_file = os.path.join(organism_input_dir, "annotations.goa")
        fasta_file = os.path.join(organism_input_dir, "proteome.fasta")
        go_filtered_file = os.path.join(organism_cache_dir, "filtered_by_go.fasta")
        mts_output_file = os.path.join(organism_cache_dir, "mts_annotations.txt")
        mts_filtered_file = os.path.join(organism_output_dir, "filtered_by_mts.fasta")

        # Step 1: Parse GO annotations
        go_annotations = parse_go_annotations(annotation_file)
        log_message(f"GO annotations parsed for {organism}.")

        # Step 2: Load proteome
        proteome = load_proteome(fasta_file)
        log_message(f"Proteome loaded for {organism}.")

        # Step 3: Filter by GO term
        filtered_by_go = filter_proteins_by_go(proteome, go_annotations, target_go_term)
        save_fasta(filtered_by_go, go_filtered_file)
        log_message(f"Proteins filtered by GO term for {organism}.")

        # Step 4: Run MitoFates
        run_perl_script(perl_script_path, go_filtered_file, "fungi", mts_output_file)
        log_message(f"MitoFates script executed for {organism}.")

        # Step 5: Filter by MTS-cleavable annotations
        mts_annotations = parse_mts_annotations(mts_output_file)
        filtered_by_mts = filter_proteins_by_mts(filtered_by_go, mts_annotations)
        save_fasta(filtered_by_mts, mts_filtered_file)
        log_message(f"Proteins filtered by MTS-cleavable annotations for {organism}.")


if __name__ == "__main__":
    main()