import argparse
import subprocess

def run_pipeline(input_file):
    # Example command to run a script with the input file
    command = ["python", "pipeline_as_command_line_argument/filter_by_go.py ", "--input", input_file]
    try:
        result = subprocess.run(command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        print("Pipeline executed successfully.")
        print(result.stdout.decode())
    except subprocess.CalledProcessError as e:
        print("Error executing pipeline:")
        print(e.stderr.decode())





if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process two command-line arguments.")
    parser.add_argument("--input", type=str, required= True, help="Input file path")
    args = parser.parse_args()
    run_pipeline(args.input)
