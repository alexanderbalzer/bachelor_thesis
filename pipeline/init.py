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
        required_keys = ["input_dir", "output_dir", "cache_dir", "target_GO_term", "mitofates_path", "cleavable"]
        for section in required_sections:
            if section not in config:
                print(f"Error: Section '{section}' not found in config.ini.")
                return False
        for key in required_keys:
            if key not in config[section]:
                print(f"Error: Key '{key}' not found in section '{section}' in config.ini.")
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


