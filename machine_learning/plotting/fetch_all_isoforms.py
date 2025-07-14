import requests
import matplotlib.pyplot as plt
import subprocess
import os
import csv
import matplotlib.patches as mpatches

# Genes of interest
genes = ["NACA", "NACAD", "NACA2"]
server = "https://rest.ensembl.org"
headers = {"Content-Type": "application/json"}

def fetch_transcripts(gene_name):
    ext = f"/lookup/symbol/homo_sapiens/{gene_name}?expand=1"
    r = requests.get(server + ext, headers=headers)
    r.raise_for_status()
    data = r.json()
    transcripts = []
    for t in data.get("Transcript", []):
        if "Translation" in t:
            transcripts.append({
                "gene": gene_name,
                "enst": t["id"],
                "ensp": t["Translation"]["id"]
            })
    return transcripts

def fetch_protein_sequence(ensp_id):
    ext = f"/sequence/id/{ensp_id}?type=protein"
    r = requests.get(server + ext, headers={"Content-Type": "text/plain"})
    r.raise_for_status()
    return r.text.strip()


def run_interproscan(
    fasta_input,
    output_file,
    interproscan_path="interproscan.sh",
    applications="Pfam,SMART",  # you can change this to "all" or specific tools
    output_format="tsv"
):
    """
    Runs InterProScan as a subprocess.

    Args:
        fasta_input (str): Path to the input FASTA file.
        output_file (str): Path to the output .tsv file.
        interproscan_path (str): Path to the interproscan.sh script.
        applications (str): Comma-separated list of applications (e.g., "Pfam,SMART").
        output_format (str): Output format, typically "tsv", "xml", or "json".

    Returns:
        None
    """
    cmd = [
        interproscan_path,
        "-i", fasta_input,
        "-f", output_format,
        "-appl", applications,
        "-dp",  # lookup protein matches
        "-goterms",
        "-pa",  # include pathway annotations
        "-o", output_file
    ]

    print(f"Running InterProScan:\n{' '.join(cmd)}")
    try:
        subprocess.run(cmd, check=True)
        print(f"✅ InterProScan completed successfully. Output saved to {output_file}")
    except subprocess.CalledProcessError as e:
        print("❌ InterProScan failed with the following error:")
        print(e)
    

all_transcripts = []
for gene in genes:
    all_transcripts.extend(fetch_transcripts(gene))

# Generate FASTA
fasta_lines = []
for t in all_transcripts:
    seq = fetch_protein_sequence(t["ensp"])
    fasta_lines.append(f">{t['enst']}|{t['ensp']}|{t['gene']}")
    fasta_lines.append(seq)

# Save FASTA
with open("naca_paralogues_proteins.fasta", "w") as f:
    f.write("\n".join(fasta_lines))

print("Saved: naca_paralogues_proteins.fasta")

fasta_input_path = "naca_paralogues_proteins.fasta"
interpro_output_path = "naca_interpro.tsv"
'''
run_interproscan(
    fasta_input=fasta_input_path,
    output_file=interpro_output_path,
    interproscan_path="/home/abalzer/interproscan-5.75-106.0-64-bit/interproscan-5.75-106.0/interproscan.sh",  # or "/path/to/interproscan.sh"
    applications="Pfam,SMART",  # use "all" for full scan
    output_format="tsv"
)
'''


domain_data = {}

with open(interpro_output_path) as f:
    reader = csv.reader(f, delimiter='\t')
    for row in reader:
        enst_id = row[0]
        enst_id = enst_id.split('|')[2] + '|' + enst_id.split('|')[0]  # Extract ENST ID
        protein_length = int(row[2])
        start = int(row[6])
        end = int(row[7])
        domain_name = row[5]
        domain_data.setdefault(enst_id, []).append((start, end, domain_name, protein_length))

# Define colors for domains

# Plotting
fig, ax = plt.subplots(figsize=(12, len(domain_data) * 0.8))
yticks = []
ylabels = []
for i, (enst, domain_info) in enumerate(domain_data.items()):
    domains = domain_info
    protein_length = domains[0][3] if domains else 0
    y = -i
    yticks.append(y)
    ylabels.append(enst)
    
    # Full protein bar
    ax.add_patch(mpatches.Rectangle((0, y - 0.2), protein_length, 0.4, color='lightgray'))
    

    displayed_domains = set()
    # Domains
    for start, end, name, _ in domains:
        # if the name contains "Nascent polypeptide", skip it
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
        ax.add_patch(mpatches.Rectangle((start, y - 0.2), end - start, 0.4,
                                             color=color, label=name))
        if name not in displayed_domains:
            displayed_domains.add(name)
            
        else:
            continue
        
        
ax.set_yticks(yticks)
ax.set_yticklabels(ylabels)

# Shift y-tick labels a little higher
for label in ax.get_yticklabels():
    label.set_y(label.get_position()[1] + 0.2)  # Adjust as needed

ax.set_xlabel("Amino Acid Position")
ax.set_title("Protein Domain Architecture of NACA, NACAD, and NACA2 Isoforms")
# Create a legend
handles, labels = ax.get_legend_handles_labels()
unique_labels = set(labels)
handles = [handles[labels.index(label)] for label in unique_labels]
ax.legend(handles, unique_labels, loc='upper right', title='Domains', fontsize='small')
ax.set_xlim(0, 2225)
ax.set_ylim(-len(domain_data) + 0.2, 1)
plt.grid(True, axis='x', linestyle='--', alpha=0.4)
plt.tight_layout()
output_path = '/home/abalzer/Documents/github_clone/bachelor_thesis/pipeline/output/output_20250617_183139_latest_ML/rna_isoforms/naca_paralogues_domains.png'
plt.savefig(output_path, dpi=300, bbox_inches='tight')
plt.show()