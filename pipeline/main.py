import init
import os
import go_filter 
import mts_filter
import heatmap
import logging
import phylogenetic_tree
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
    delete_cache = config['DEFAULT']['delete_cache']
    save_filtered_proteins = config['DEFAULT'].getboolean('save_filtered_proteins', True)
    create_heatmap = config['DEFAULT'].getboolean('create_heatmap', True)
    heatmap_type = config['DEFAULT'].get('heatmap_type', 'hgt')
    create_phylogenetic_tree = config['DEFAULT'].getboolean('create_phylogenetic_tree', True)
    phylo_tree_type = config['DEFAULT'].get('phylo_tree_type', 'absolute')
    phylo_tree_method = config['DEFAULT'].get('phylo_tree_method', 'pearson')
    phylo_tree_algorithm = config['DEFAULT'].get('phylo_tree_algorithm', 'upgma')
    save_newick = config['DEFAULT'].getboolean('save_newick', True)


    # Log the start of the pipeline
    log_message("Starting the pipeline...")

    # Read organism names from FASTA file names
    fasta_files = [f for f in os.listdir(input_dir) if f.endswith(".fasta")]
    organism_names = [os.path.splitext(f)[0] for f in fasta_files]
    log_message(f"Organism names extracted: {', '.join(organism_names)}")
    
    # Filter proteins by GO term
    go_filter.run(organism_names, input_dir, cache_dir, target_go_term)
    log_message("GO term filtering completed.")

    #read the flaglist for MitoFates as a dictionairy
    flaglist = {}
    with open(os.path.join(input_dir, "flaglist.txt"), "r") as file:
        for line in file:
            key, value = line.strip().split(":")
            flaglist[key.strip()] = value.strip()
    log_message(f"Flaglist loaded: {flaglist}")

    # Filter proteins by MTS-cleavable probability
    mts_filter.run(organism_names, cache_dir, output_dir, cleavable, perl_script_path, flaglist, delete_cache, threshold)
    log_message("MitoFates filtering completed.")

    if create_heatmap or create_phylogenetic_tree:
        heatmap.run(organism_names, cache_dir, output_dir, create_heatmap, heatmap_type, create_phylogenetic_tree, phylo_tree_type)
        if create_heatmap:
            log_message("Heatmap creation completed.")

    if create_phylogenetic_tree:
        phylogenetic_tree.run(organism_names, cache_dir, output_dir, phylo_tree_method, phylo_tree_algorithm, save_newick)
        log_message("Phylogenetic tree creation completed.")

    if not save_filtered_proteins:
        for organism in organism_names:
            os.remove(os.path.join(output_dir, f"{organism}_filtered_by_go_and_mts.fasta"))
            log_message(f"Deleted filtered proteins file for {organism}")
    else:
        log_message("Filtered proteins files are retained.")


    log_message("Pipeline completed successfully.")


if __name__ == "__main__":
    main()