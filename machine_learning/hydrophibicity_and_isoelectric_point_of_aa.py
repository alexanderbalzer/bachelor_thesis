import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from adjustText import adjust_text

# Data for amino acids
amino_acids = [
    'A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I',
    'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'
]
hydrophobicity = [
    1.8, -4.5, -3.5, -3.5, 2.5, -3.5, -3.5, -0.4, -3.2, 4.5,
    3.8, -3.9, 1.9, 2.8, -1.6, -0.8, -0.7, -0.9, -1.3, 4.2
]
isoelectric_points = [
    6.01, 10.76, 5.41, 2.77, 5.07, 3.22, 5.65, 5.97, 7.59, 6.02,
    5.98, 9.74, 5.74, 5.48, 6.30, 5.68, 5.60, 5.89, 5.66, 5.96
]
NatAD = ["A", "C", "T", "S", "V", "G"]
NatB = ["D", "E", "N", "Q"]
NatCE = ["L", "I", "F", "Y", "K"]

polar_aa = ['S', 'T', 'N', 'H', 'Q', 'G']
speci_aa = ['P', 'C']
apolar_aa = ['A', 'L', 'V', 'I', 'M']
charged_aa = ['E', 'D', 'K', 'R']
aromatic_aa = ['W', 'Y', 'F']

# Assign colors and markers
colors = []
markers = []
for aa in amino_acids:
    if aa in NatAD:
        colors.append('tab:blue')
    elif aa in NatB:
        colors.append('tab:orange')
    elif aa in NatCE:
        colors.append('tab:green')
    else:
        colors.append('tab:gray')

    # Assign marker by property
    if aa in polar_aa:
        markers.append('o')      # circle
    elif aa in speci_aa:
        markers.append('s')      # square
    elif aa in apolar_aa:
        markers.append('^')      # triangle_up
    elif aa in charged_aa:
        markers.append('D')      # diamond
    elif aa in aromatic_aa:
        markers.append('*')      # star
    else:
        markers.append('x')      # cross

# Create the scatter plot
plt.figure(figsize=(10, 6))
texts = []
for i, aa in enumerate(amino_acids):
    plt.scatter(hydrophobicity[i], isoelectric_points[i], 
                c=colors[i], marker=markers[i], s=120, alpha=0.7, edgecolor='black')
    texts.append(
        plt.text(hydrophobicity[i], isoelectric_points[i], aa, fontsize=15, ha='right', color=colors[i])
    )

# Adjust text to avoid overlaps
adjust_text(texts, arrowprops=dict(arrowstyle='-', color='gray', lw=0.5))

plt.xlabel('Hydrophobicity')
plt.ylabel('Isoelectric Point')
plt.title('Isoelectric Points vs Hydrophobicity of Amino Acids')
plt.grid(True)

# Custom legend for color (Nat types)
legend_elements = [
    Line2D([0], [0], marker='o', color='w', label='NatAD', markerfacecolor='tab:blue', markersize=10),
    Line2D([0], [0], marker='o', color='w', label='NatB', markerfacecolor='tab:orange', markersize=10),
    Line2D([0], [0], marker='o', color='w', label='NatCE', markerfacecolor='tab:green', markersize=10),
    Line2D([0], [0], marker='o', color='w', label='Other', markerfacecolor='tab:gray', markersize=10)
]
plt.legend(handles=legend_elements, title="N-terminal Acetylation Type", loc='upper left', bbox_to_anchor=(1, 1))

# Custom legend for marker shape (property)
marker_legend = [
    Line2D([0], [0], marker='o', color='k', label='Polar', markerfacecolor='w', markersize=10),
    Line2D([0], [0], marker='s', color='k', label='Special', markerfacecolor='w', markersize=10),
    Line2D([0], [0], marker='^', color='k', label='Apolar', markerfacecolor='w', markersize=10),
    Line2D([0], [0], marker='D', color='k', label='Charged', markerfacecolor='w', markersize=10),
    Line2D([0], [0], marker='*', color='k', label='Aromatic', markerfacecolor='w', markersize=10),
    Line2D([0], [0], marker='x', color='k', label='Other', markerfacecolor='w', markersize=10)
]
plt.legend(handles=legend_elements + marker_legend, title="Legend", loc='upper left', bbox_to_anchor=(1, 1))

plt.tight_layout()
plt.savefig("pipeline/output/output_20250519_142700_machine_learning_human/Homo_sapiens/hydrophobicity_vs_isoelectric_point.png", dpi=300)

plt.show()