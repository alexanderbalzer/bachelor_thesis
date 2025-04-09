import init
import os
import go_filter 
import mts_filter
import logging
from utils import log_message

def main():
    """
    Main function to run the pipeline.
    """

    logging.basicConfig(level=logging.INFO)
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
    perl_script_path = os.path.abspath(config['DEFAULT']['mitofates_path'])
    cleavable = config['DEFAULT'].get('cleavable', 'No')
    threshold = config['DEFAULT'].getfloat('threshold', 0.9)
    save_filtered_proteins = config['DEFAULT'].getboolean('save_filtered_proteins', False)
    save_hgt_array = config['DEFAULT'].getboolean('save_hgt_array', False)
    delete_cache = config['DEFAULT'].getboolean('delete_cache', False)


    # Log the start of the pipeline
    log_message("Starting the pipeline...")

    # Read organism names from FASTA file names
    fasta_files = [f for f in os.listdir(input_dir) if f.endswith(".fasta")]
    organism_names = [os.path.splitext(f)[0] for f in fasta_files]
    log_message(f"Organism names extracted: {', '.join(organism_names)}")
    
    # Filter proteins by GO term
    go_filter.run(organism_names, input_dir, cache_dir, output_dir, target_go_term)
    log_message("GO term filtering completed.")

    #read the flaglist for MitoFates as a dictionairy
    flaglist = {}
    with open(os.path.join(input_dir, "flaglist.txt"), "r") as file:
        for line in file:
            parts = line.strip().split("\t")
            if len(parts) == 2:
                flaglist[parts[0]] = parts[1]
    log_message(f"Flaglist loaded: {flaglist}")

    # Filter proteins by MTS-cleavable probability






if __name__ == "__main__":
    main()