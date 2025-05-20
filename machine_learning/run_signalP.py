import subprocess
import os

def run_signalp(fasta_file, output_dir, signalp_path="signalp6"):
    """
    Run SignalP 6 on a given FASTA file.

    Args:
        fasta_file (str): Path to the input FASTA file.
        output_dir (str): Directory to save the SignalP output.
        signalp_path (str): Path to the SignalP 6 executable (default: "signalp6").
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    output_file = os.path.join(output_dir, "signalp_output.txt")

    command = [
        signalp_path,
        "--fastafile", fasta_file,
        "--organism", "eukarya",
        "--output_dir", output_dir,
        "--format", "txt",
        "--mode", "fast"
    ]

    try:
        subprocess.run(command, check=True)
        print(f"SignalP 6 analysis completed. Results saved to {output_file}")
    except subprocess.CalledProcessError as e:
        print(f"Error running SignalP 6: {e}")
    except FileNotFoundError:
        print(f"SignalP 6 executable not found at {signalp_path}. Please check the installation.")

# Example usage
if __name__ == "__main__":
    fasta_file = "/home/abalzer/Documents/github_clone/bachelor_thesis/pipeline/input/Homo_sapiens.fasta" 
    output_dir = "/home/abalzer/Documents/github_clone/bachelor_thesis/pipeline/output/output_20250519_142700_machine_learning_human/Homo_sapiens" 
    run_signalp(fasta_file, output_dir)