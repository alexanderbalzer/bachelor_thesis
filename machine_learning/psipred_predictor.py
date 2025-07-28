import os
import subprocess
import tempfile
from pathlib import Path

psipred_path = os.getenv("PSIPRED_PATH", "/home/abalzer/psipred")
runpsipred_single = os.path.join(psipred_path, "runpsipred_single")

def run_psipred(sequence: str) -> float:
    """Run PSIPRED on a peptide sequence and return % helices."""

    # Handle ACE- acetylation prefix
    sequence = sequence.replace("ACE-", "").upper()
    if not sequence.isalpha():
        raise ValueError("Invalid sequence format.")

    # Set up temporary directory
    with tempfile.TemporaryDirectory() as tmpdir:
        fasta_path = Path(tmpdir) / "input.fasta"
        with open(fasta_path, "w") as f:
            f.write(">query\n" + sequence + "\n")

        # Run psipred pipeline (must be configured with proper paths)
        cmd = f"{runpsipred_single} {fasta_path}"
        result = subprocess.run(cmd, shell=True, cwd=tmpdir, capture_output=True)

        if result.returncode != 0:
            raise RuntimeError("PSIPRED failed:\n" + result.stderr.decode())

        horiz_path = Path(tmpdir) / "input.horiz"
        if not horiz_path.exists():
            raise FileNotFoundError("PSIPRED output .horiz file not found.")

        # Parse .horiz output to get secondary structure prediction
        helix_count = 0
        total = 0
        with open(horiz_path) as f:
            for line in f:
                if line.startswith("Pred:"):
                    pred_line = line.strip().split()[1]
                    helix_count += pred_line.count('H')
                    total += len(pred_line)

        if total == 0:
            raise ValueError("No residues processed.")
        return round(helix_count / total * 100, 2)

if __name__ == "__main__":
    # Example usage
    sequence = "ACE-MKLLLILLLLILLLL"
    try:
        helix_percentage = run_psipred(sequence)
        print(f"Helix percentage: {helix_percentage}%")
    except Exception as e:
        print(f"Error: {e}")