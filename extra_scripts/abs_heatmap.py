import numpy as np
from collections import Counter
from Bio import SeqIO
import matplotlib.pyplot as plt


def classify_natc_substrate(second_as):
    if second_as in ["A", "C", "T", "S", "V", "G"]:
        return "NatA/D"
    elif second_as in ["D", "E", "N", "Q"]:
        return "NatB"
    elif second_as in ["L", "I", "F", "Y", "K"]:
        return "NatC/E"
    else:
        return "Other"
    


# Path to the FASTA file
fasta_file = "/home/abalzer/Documents/github_clone/bachelor_thesis/pipeline/output/output_20250506_111626/human/human_filtered_by_GO_cleavable_mts.fasta"

# Initialize a counter for amino acids
amino_acid_counts = Counter()

# Parse the FASTA file and count amino acids for each position
position_counts = [Counter() for _ in range(20)]  # Initialize counters for each position

for record in SeqIO.parse(fasta_file, "fasta"):
    for i, amino_acid in enumerate(record.seq[:20]):  # Only consider the first 20 positions
        position_counts[i][amino_acid] += 1

# Get the amino acids (alphabetically sorted) and their frequencies for each position
amino_acids = sorted(set(aa for pos in position_counts for aa in pos.keys()))
heatmap_data = np.zeros((len(amino_acids), 20))

for i, pos_count in enumerate(position_counts):
    for j, aa in enumerate(amino_acids):
        heatmap_data[j, i] = pos_count[aa]

# Reorder the heatmap data to match the specified amino acid order
amino_acid_order = ["D", "E", "N", "Q", "Y", "H", "K", "R", "M", "L", "F", "I", "W", "S", "A", "T", "C", "P", "G", "V"]
reordered_indices = [amino_acids.index(aa) for aa in amino_acid_order]
heatmap_data = heatmap_data[reordered_indices, :]
amino_acids = amino_acid_order

# Create a heatmap
fig, ax = plt.subplots(figsize=(10, 8))
heatmap_data = heatmap_data[:, 1:]  # Exclude the first position
heatmap_data = np.log10(heatmap_data + 1)  # Log scale for better visualization
cax = ax.imshow(heatmap_data, cmap="viridis", aspect="auto")

# Add labels and colorbar
ax.set_yticks(range(len(amino_acids)))
ax.set_yticklabels(amino_acids)
ax.set_xticks(range(19))
ax.set_xticklabels(range(1, 20))  # Positions 1 to 20
plt.colorbar(cax, orientation="vertical", label="Frequency")

# Title and display
plt.title("Heatmap of Amino Acid Frequencies by Position (First 20 Positions)")
plt.xlabel("Position")
plt.ylabel("Amino Acids")
plt.tight_layout()
plt.show()