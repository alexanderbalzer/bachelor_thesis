import requests
import matplotlib.pyplot as plt
import subprocess
import os
import csv
import matplotlib.patches as mpatches
import pandas as pd


"""
This script uses InterPro domain annotations as input and plots the domain architecture of the isoforms in a horizontal bar chart.
It also filters the isoforms based on a predefined list and visualizes the domain architecture.
Input is a TSV file with InterPro annotations and a list of isoforms.
"""

# Genes of interest
genes = ["NACA", "NACAD", "NACA2", "BTF3", "BTF3L4"]
server = "https://rest.ensembl.org"
headers = {"Content-Type": "application/json"}

    

all_transcripts = []
working_dir = os.getcwd()  # Use current working directory
interpro_output_path = os.path.join(working_dir, "Isoform plotting", "naca_interpro.tsv")

domain_data = {}

with open(interpro_output_path) as f:
    reader = csv.reader(f, delimiter='\t')
    for row in reader:
        enst_id = row[0]
        enst_id = enst_id.split('|')[0]  # Extract ENST ID
        protein_length = int(row[2])
        start = int(row[6])
        end = int(row[7])
        domain_name = row[5]
        domain_data.setdefault(enst_id, []).append((start, end, domain_name, protein_length))

isoform_list = pd.read_csv(os.path.join(working_dir, "Isoform plotting", "isoform_list.txt"), header=None).squeeze().tolist()

# use only isoforms that are in the isoform_list
domain_data = {enst: domains for enst, domains in domain_data.items() if enst in isoform_list}
# order the domain_data by the isoform_list
domain_data = {enst: domain_data[enst] for enst in isoform_list if enst in domain_data}

# Define colors for domains

# Plotting
fig, ax = plt.subplots(figsize=(len(domain_data) * 0.8, 9))  # Swap figsize

xticks = []
xlabels = []
for i, (enst, domain_info) in enumerate(domain_data.items()):
    domains = domain_info
    protein_length = domains[0][3] if domains else 0
    x = i
    xticks.append(x)
    xlabels.append(enst)
    
    # Full protein bar (now vertical)
    ax.add_patch(mpatches.Rectangle((x - 0.2, 0), 0.4, protein_length, color='lightgray'))
    
    displayed_domains = set()
    # Domains (now vertical)
    for start, end, name, _ in domains:
        if "disorder" in name.lower():
            name = "Disordered"
            color = 'orange'
        elif "nascent" in name.lower():
            continue
        elif "uba" in name.lower():
            name = "UBA"
            color = 'lightgreen'
        elif "nac" in name.lower():
            name = "NAC"
            color = 'darkgreen'
        else:
            continue
        ax.add_patch(mpatches.Rectangle((x - 0.2, start), 0.4, end - start, color=color, label=name))
        if name not in displayed_domains:
            displayed_domains.add(name)
        else:
            continue

# Collect legend handles for unique domains
legend_dict = {
    "Disordered": 'orange',
    "UBA": 'lightgreen',
    "NAC": 'darkgreen',
    "Unannotated ": 'lightgray'
}
handles = [mpatches.Patch(color=color, label=name) for name, color in legend_dict.items()]
ax.legend(handles=handles, loc='lower right', title='Domain', fontsize='x-large', title_fontsize='x-large')

ax.set_xticks(xticks)
xlabels = "" * len(xticks)  # Initialize xlabels as empty strings
ax.set_xticklabels(xlabels, rotation=45, ha='right')
ax.set_xlabel("")
ax.xaxis.set_ticks_position('top')
ax.xaxis.set_label_position('top')

ax.set_ylabel("Amino Acid Position")
ax.set_title("")

ax.set_xlim(-0.5, len(domain_data) - 0.5)
ax.set_ylim(0, 2100)
ax.invert_yaxis()  # Mirror the plot on the x-axis

plt.grid(True, axis='y', linestyle='--', alpha=0.4)
plt.tight_layout()
output_path = os.path.join(working_dir, "domain_architecture_flipped_mirrored.png")
plt.savefig(output_path, dpi=300, bbox_inches='tight')