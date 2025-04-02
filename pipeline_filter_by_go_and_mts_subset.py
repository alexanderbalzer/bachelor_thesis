import subprocess
import sys



name = ["human"]

scripts = [
    "filter_mitochondrial_proteins_from_proteome.py",
    "subset_MTS_targetting_sequence.py",
    "Heatmap_of_first_AS_of_Proteins.py"
]
for script in scripts:
    print(f"Running {script} for {name}")
    subprocess.run(["python", script], check=True)
    print(f"Finished running {script} for {name}")









