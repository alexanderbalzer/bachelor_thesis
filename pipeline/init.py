import configparser
import os

def check_config():
    """
    Check if the config.ini file exists.
    """
    config_path = os.path.join(os.path.dirname(__file__), "config.ini")
    if not os.path.exists(config_path):
        print("Error: config.ini file not found.")
        return False
    else:
    # Check for required sections and keys
        config = configparser.ConfigParser()
        config.read(config_path)
        section = "DEFAULT"
        required_sections = ["DEFAULT"]
        required_keys = ["input_dir", "output_dir", "cache_dir", "target_GO_term", 
                         "mitofates_path", "cleavable", "threshold", "save_filtered_proteins", 
                         "delete_cache", "create_heatmap", "create_phylogenetic_tree", 
                         "type", "phylo_tree_method", "phylo_tree_algorithm", "save_newick"]
        # Check if the required sections and keys are present
        for section in required_sections:
            if section not in config:
                print(f"Error: Section '{section}' not found in config.ini.")
                return False
        for key in required_keys:
            if key not in config[section]:
                print(f"Error: Key '{key}' not found in section '{section}' in config.ini.")
                return False
        
        # check if the perl script path exists
        perl_script_path = config[section]['mitofates_path']
        if not os.path.isabs(perl_script_path) or not os.path.exists(perl_script_path):
            print(f"Error: Perl script '{perl_script_path}' not found or the path is not absolute.")
            return False
        #check if the output_dir and cache_dir exist and if not create them
        if not os.path.exists(config[section]['output_dir']):
            os.makedirs(config[section]['output_dir'])
        if not os.path.exists(config[section]['cache_dir']):
            os.makedirs(config[section]['cache_dir'])

        return True
    

def load_config():
    """
    Load the configuration from config.ini.
    """
    config = configparser.ConfigParser()
    config.read(os.path.join(os.path.dirname(__file__), "config.ini"))
    return config

def transform_labels_to_names(organism_names):
    """
    Transform organism labels to names 
    """
    names = []
    for label in organism_names:
        # Replace underscores with spaces
        name = label.replace("_", " ")
        # Capitalize the first letter of each word
        name = name.title()
        # Add the transformed name to the list
        names.append(name)

    return names
