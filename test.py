from Bio import SeqIO
import pandas as pd
import pandas as pd
from Bio import SeqIO
import os
import logging
import subprocess

logging.basicConfig(level=logging.INFO)


def run_perl_script(perl_script_path, input_file, flag, output_file):
    with open(output_file, "w") as output_file:
        subprocess.run(
            ["perl", perl_script_path, input_file, flag],
            stdout=output_file,
            stderr=subprocess.PIPE,
            check=True
    )

def check_file_exists(file_path):
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"File not found: {file_path}")



name = "human"
perl_script_path = "/home/abalzer/Downloads/MitoFates_1.2/MitoFates/MitoFates.pl"
#input_file = "pipeline/cache/filtered_proteins_by_GO_for_" + str(name) + ".fasta"
#output_file = "pipeline/cache/mito_fates_for_" + str(name) + ".cgi"
input_file = "pipeline/cache/example.fasta"
output_file = "pipeline/cache/example.cgi"



run_perl_script(perl_script_path, input_file, "fungi", output_file)


logging.info(f"Perl script executed successfully. Output saved to {output_file}")







