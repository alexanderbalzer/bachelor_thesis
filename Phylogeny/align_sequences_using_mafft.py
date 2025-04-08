import subprocess

def align_sequences(input_fasta, output_fasta):
    try:
        # Run MAFFT command
        result = subprocess.run(
            ["mafft", "--auto", input_fasta],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        
        # Check for errors
        if result.returncode != 0:
            print("Error running MAFFT:")
            print(result.stderr)
            return
        
        # Write the aligned sequences to the output file
        with open(output_fasta, "w") as output_file:
            output_file.write(result.stdout)
        
        print(f"Alignment completed. Results saved to {output_fasta}")
    except FileNotFoundError:
        print("MAFFT is not installed")

# Example usage
input_fasta = "Phylogeny/input/18s.fasta"  # Replace with your input FASTA file
output_fasta = "Phylogeny/output/aligned_sequences.fasta"  # Replace with your desired output file
align_sequences(input_fasta, output_fasta)