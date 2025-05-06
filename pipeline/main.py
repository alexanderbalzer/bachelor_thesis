import init
import os
import go_filter 
import mts_filter
import heatmap
import logging
import phylogenetic_tree
import Nat_GOterm_link
import logoplot
import Heatmap_one_organism_absolute_or_HGT
from utils import log_message
from datetime import datetime
import pandas as pd
import shutil

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
    delete_cache = config['DEFAULT']['delete_cache']
    save_filtered_proteins = config['DEFAULT'].getboolean('save_filtered_proteins', True)
    create_heatmap = config['DEFAULT'].getboolean('create_heatmap', True)
    heatmap_type = config['DEFAULT'].get('heatmap_type', 'hgt')
    create_phylogenetic_tree = config['DEFAULT'].getboolean('create_phylogenetic_tree', True)
    type = config['DEFAULT'].get('type', 'absolute')
    phylo_tree_method = config['DEFAULT'].get('phylo_tree_method', 'pearson')
    phylo_tree_algorithm = config['DEFAULT'].get('phylo_tree_algorithm', 'upgma')
    save_newick = config['DEFAULT'].getboolean('save_newick', True)
    run_from_scratch = config['DEFAULT'].getboolean('run_from_scratch', False)
    reference = config['DEFAULT'].get('reference', 'subset')
    create_logoplot = config['DEFAULT'].getboolean('create_logoplot', True)

    log_message("Configuration loaded successfully.")


    # Log the start of the pipeline
    log_message("Starting the pipeline...")

    # Read organism names from FASTA file names
    fasta_files = [f for f in os.listdir(input_dir) if f.endswith(".fasta")]
    organism_names = [os.path.splitext(f)[0] for f in fasta_files]
    log_message(f"Organism names extracted: {', '.join(organism_names)}")

    # Read the last run timestamp from the file cache_of_last_run.txt
    last_run_file = os.path.join(os.getcwd(), "cache_of_last_run.txt")
    if os.path.exists(last_run_file):
        with open(last_run_file, "r") as file:
            last_run = file.read().strip()
        log_message(f"Last run timestamp loaded: {last_run}")
    else:
        last_run = None
        run_from_scratch = True
        log_message("No previous run timestamp found.")

    # Create a folder with a timestamp for the cache and the output files
    now = datetime.now()
    timestamp = now.strftime("%Y%m%d_%H%M%S")
    cache_dir = os.path.join(cache_dir, f"cache_{timestamp}/")
    with open(last_run_file, "w") as file:
        file.write(cache_dir)
    output_dir = os.path.join(output_dir, f"output_{timestamp}/")
    os.makedirs(cache_dir, exist_ok=True)
    log_message(f"Cache directory created: {cache_dir}")
    os.makedirs(output_dir, exist_ok=True)
    log_message(f"Output directory created: {output_dir}")
    # copy the config file to the output directory
    config_file_path = "pipeline/config.ini"
    if os.path.exists(config_file_path):
        shutil.copy(config_file_path, output_dir)
        log_message(f"Config file copied to output directory: {output_dir}")

    # Initialize a DataFrame to keep track of the amount of proteins per step
    amount_of_proteins_per_step = pd.DataFrame(index=["Start", "Mitochondrial", "Mitochondrial with MTS"], columns=organism_names)

    # Filter proteins by GO term
    amount_of_proteins_per_step = go_filter.run(organism_names, input_dir, output_dir, target_go_term, run_from_scratch, amount_of_proteins_per_step, last_run)
    log_message("GO term filtering completed.")

    #read the flaglist for MitoFates as a dictionairy
    flaglist = {}
    flaglist_path = os.path.join(os.path.dirname(__file__), "flaglist.txt")
    if not os.path.exists(flaglist_path):
        print("Error: flaglist.txt file not found.")
        return False
    with open(flaglist_path, "r") as file:
        for line in file:
            key, value = line.strip().split(":")
            flaglist[key.strip()] = value.strip()
    log_message(f"Flaglist loaded: {flaglist}")

    # Filter proteins by MTS-cleavable probability
    amount_of_proteins_per_step = mts_filter.run(organism_names, output_dir, cleavable, perl_script_path, flaglist, delete_cache, threshold, run_from_scratch, amount_of_proteins_per_step, last_run)
    log_message("MitoFates filtering completed.")
    # save the amount of proteins per step to a file
    amount_of_proteins_per_step.to_csv(os.path.join(output_dir, "amount_of_proteins_per_step.csv"), index=False)
    log_message("Amount of proteins per step saved.")

    if create_heatmap or create_phylogenetic_tree:
        heatmap.run(organism_names, input_dir, cache_dir, output_dir, create_heatmap, heatmap_type, create_phylogenetic_tree, type, reference)
        if create_heatmap:
            log_message("Heatmap creation completed.")

    if create_phylogenetic_tree:
        phylogenetic_tree.run(organism_names, cache_dir, output_dir, phylo_tree_method, phylo_tree_algorithm, save_newick)
        log_message("Phylogenetic tree creation completed.")
    
    log_message("Creating heatmaps for each organism.")
    Heatmap_one_organism_absolute_or_HGT.run(organism_names, input_dir, output_dir, heatmap_type)
    log_message("Heatmap creation for each organism completed.")
    
    log_message("Creating Nat GO term link.")
    Nat_GOterm_link.run(organism_names, input_dir, output_dir)
    log_message("Nat GO term link creation completed.")
    
    log_message("creating logoplot")
    logoplot.run(organism_names, output_dir)
    log_message("Logoplot creation completed.")



    if not save_filtered_proteins:
        for organism in organism_names:
            os.remove(os.path.join(output_dir, f"{organism}_filtered_by_go_and_mts.fasta"))
            log_message(f"Deleted filtered proteins file for {organism}")
    else:
        log_message("Filtered proteins files are retained.")
    
    if delete_cache.lower() == "yes" and not save_filtered_proteins:
        os.rmdir(cache_dir)
        log_message(f"Deleted cache directory: {cache_dir}")

    # Log the completion of the pipeline
    log_message("Pipeline completed successfully.")




if __name__ == "__main__":
    main()