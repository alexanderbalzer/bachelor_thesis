import pandas as pd
import seaborn as sns
import numpy as np
import os
import matplotlib.pyplot as plt

# === 1. Eingabedatei ===
working_dir = 'pipeline/output/output_20250617_183139_latest_ML/RNA expression consensus'
input_file = os.path.join(working_dir, 'rna_tissue_consensus.tsv')

# === 2. Ziel-Gene definieren ===
target_genes = ["NACA", "NACA2", "NACAD", 
                "BTF3", "METAP1", "METAP2", 
                "HYPK", "NAA10", "NAA50", "NAA15",
                "NAA20", "NAA25", 
                "NAA38", "NAA30", "NAA35"]
target_genes = ["NACA", "NACA2", "NACAD", 'BTF3', "BTF3L4"]

# === 3. Datei einlesen ===
df = pd.read_csv(input_file, sep="\t")

# === 4. Nur Ziel-Gene filtern ===
filtered_df = df[df["Gene name"].isin(target_genes)]
filtered_df = filtered_df[["Gene name", "Tissue", "nTPM"]]


# === 5. Pivot-Tabelle erstellen ===
expression_matrix = filtered_df.pivot_table(index="Gene name", columns="Tissue", values="nTPM", )

# === 6. Optional: fehlende Werte mit 0 oder NaN auffüllen ===
expression_matrix.fillna(0, inplace=True)  # oder: dropna(axis=1) für nur vollständige Gewebe

# === 6b. Zeilenweise Normalisierung (Maximum jeder Zeile ist 1) ===
#expression_matrix = expression_matrix.div(expression_matrix.max(axis=1), axis=0)
#expression_matrix = expression_matrix.div(expression_matrix.sum(axis=0), axis=1)
#expression_matrix = expression_matrix.div(expression_matrix.sum(axis=1), axis=0)
#expression_matrix = expression_matrix.div(expression_matrix.mean(axis=1), axis=0)
expression_matrix = np.log10(expression_matrix + 1)

# Add the highest (absolute, not normalized) expression value to the y-axis labels
max_expression = filtered_df.pivot_table(index="Gene name", values="nTPM", aggfunc="max")
# Format: "GENE (max: 12.3)"
'''new_index = [
    f"{gene}                   "
    for gene in expression_matrix.index
]
expression_matrix.index = new_index'''

# === 7. Clustermap zeichnen ===

g = sns.clustermap(
    expression_matrix.T,
    cmap="PuBuGn",
    figsize=(6, 8),
    metric="euclidean",
    yticklabels=True,  # Ensure all x-tick labels are shown
    #vmin=0, vmax=0.1,
    cbar_pos=(0.7, 0.85, 0.2, 0.03), # Position of the colorbar (x, y, width, height)
    cbar_kws={"label": "log10(nTPM)", 'orientation': "horizontal", 'ticks': np.arange(0, 4, 1)},
    row_cluster=True,  # Cluster rows (genes)
    col_cluster=False,  # Cluster columns (tissues)
)
plt.setp(g.ax_heatmap.get_xticklabels(), rotation=45, ha='right')  # Rotate and right-align
plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0, ha='left')  # Rotate y-tick labels to horizontal
g.ax_heatmap.set_ylabel("")

# Define your tissue groups
brain = {"cerebral cortex", "cerebellum", "basal ganglia", "hippocampal formation", "hypothalamus", "midbrain", "amygdala", "choroid plexus", "spinal cord", "retina"}
glands = {"thyroid gland", "pituitary gland", "adrenal gland", "parathyroid gland"}

# Color x-axis labels
xticklabels = g.ax_heatmap.get_yticklabels()
for label in xticklabels:
    tissue = label.get_text()
    if tissue in brain:
        label.set_color('darkgoldenrod')
    elif tissue in glands:
        label.set_color('purple')
    elif tissue == "lung":
        label.set_color('forestgreen')
    elif tissue in {"salivary gland", "esophagus", "tongue"}:
        label.set_color('orange')
    elif tissue in {"duodenum", "small intestine", "colon", "rectum"}:
        label.set_color('blue')
    elif tissue in {"liver", "gallbladder"}:
        label.set_color('orchid')
    elif tissue == "pancreas":
        label.set_color('green')
    elif tissue in {"testis", "epididymis", "prostate", "seminal vesicle"}:
        label.set_color('deepskyblue')
    elif tissue in {"vagina", "ovary," "cervix", "endometrium", "placenta", "fallopian tube", "breast"}:
        label.set_color('pink')
    elif tissue in {"heart muscle", "smooth muscle", "skeletal muscle"}:
        label.set_color('brown')
    elif tissue in {"adipose tissue"}:
        label.set_color('aqua')
    elif tissue in {"skin"}:
        label.set_color('tan')
    elif tissue in {"appendix", "spleen", "lymph node", "thymus", "tonsil", "bone marrow"}:
        label.set_color('red')
    else:
        label.set_color('black')

yticklabels = g.ax_heatmap.get_xticklabels()
for label in yticklabels:
    protein = label.get_text().split(" (")[0]  # Extract gene name before parentheses
    if protein in {"HYPK", "NAA15"}:
        label.set_color('green')
    if protein == "NAA10":
        label.set_color('lightgreen')
    if protein == "NAA50":
        label.set_color('darkgreen')
    if protein in {"NAA20", "NAA25"}:
        label.set_color('blue')
    if protein in {"NAA30", "NAA35", "NAA38"}:
        label.set_color('darkred')
    


# Add grid lines to the heatmap
ax = g.ax_heatmap
num_y, num_x = expression_matrix.shape
for x in range(num_y + 1):
    ax.axvline(x , color='black', linewidth=0.1)
for y in range(num_x + 1):
    ax.axhline(y, color='black', linewidth=0.1)

plt.savefig(os.path.join(working_dir, "clustermap_ntpm.png"), dpi=300)

# === 8. Optional: Tabelle speichern ===
expression_matrix.to_csv(os.path.join(working_dir, "nat_expression_matrix.csv"))
