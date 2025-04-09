import argparse
import subprocess
import os
import configparser

def run_pipeline(input_directory):
    # Example command to run a script with the input file
    command = ["python", "pipeline_as_command_line_argument/filter_by_go.py ", "--input", input_directory]
    try:
        result = subprocess.run(command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        print("filtered_by_GO_successfully")
        print(result.stdout.decode())
    except subprocess.CalledProcessError as e:
        print("Error executing pipeline:")
        print(e.stderr.decode())





if __name__ == "__main__":
    # Check for a config file in the same directory
    config_file_path = os.path.join(os.path.dirname(__file__), "config.ini")
    config = configparser.ConfigParser()

    if os.path.exists(config_file_path):
        config.read(config_file_path)
        print("Config file loaded successfully.")
    else:
        print("No config file found. Please provide.")

    # Check if the input directory is specified in the config file
    if 'DEFAULT' in config and 'input_directory' in config['DEFAULT']:
        input_directory = config['DEFAULT']['input_directory']
    else:
        print("No input directory specified in the config file.")
        input_directory = None

    # Parse command line arguments
    parser = argparse.ArgumentParser(description="tba")
    parser.add_argument("--input", type=str, required= True, help="Input file path")
    args = parser.parse_args()
    run_pipeline(args.input)

