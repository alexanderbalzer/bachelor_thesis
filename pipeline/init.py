import os
import configparser



def check_config():
    """
    Check if the config file exists and load it.
    """
    clearance = True
    # Check for a config file in the same directory
    config_file_path = os.path.join(os.path.dirname(__file__), "config.ini")
    config = configparser.ConfigParser()

    if os.path.exists(config_file_path):
        config.read(config_file_path)
        print("Config file loaded successfully.")
    else:
        print("No config file found. Please provide.")
        clearance = False
    # Check for required sections in the config file
    required_sections = ["DEFAULT"]
    for section in required_sections:
        if section not in config:
            print(f"Missing section: {section} in config file.")
            clearance = False
    # Check for required keys in the DEFAULT section
    required_keys = [
        "input_dir",
        "mitofates_path",
        "target_go_term",
        "cleavable",
        "threshold",
        "cache_dir",
        "output_dir"
    ]
    for key in required_keys:
        if key not in config["DEFAULT"]:
            print(f"Missing key: {key} in config file.")
            clearance = False
    # Check for required directories
    required_dirs = [
        config["DEFAULT"]["input_dir"],
        config["DEFAULT"]["mitofates_path"],
    ]
    for dir_path in required_dirs:
        if not os.path.exists(dir_path):
            print(f"Directory does not exist: {dir_path}")
            clearance = False
    return clearance

def load_config():
    """
    Load the configuration file.
    """
    config_file_path = os.path.join(os.path.dirname(__file__), "config.ini")
    config = configparser.ConfigParser()
    config.read(config_file_path)
    return config


