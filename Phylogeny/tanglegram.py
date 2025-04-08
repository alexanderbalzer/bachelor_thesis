import matplotlib.pyplot as plt
from io import StringIO
from Bio import Phylo


# Read Newick strings from two files
with open("Phylogeny/output/phylogenetic_tree_from_HGT.newick", "r") as file1, open("Phylogeny/output/phylogenetic_tree_from_16s.newick", "r") as file2:
    newick1 = file1.read().strip()
    newick2 = file2.read().strip()

# Parse the trees using Biopython
tree1 = Phylo.read(StringIO(newick1), "newick")
tree2 = Phylo.read(StringIO(newick2), "newick")

# Create a mapping of leaf names to positions for both trees
def get_leaf_positions(tree):
    leaf_positions = {}
    for i, leaf in enumerate(tree.get_terminals()):
        leaf_positions[leaf.name] = i
    return leaf_positions

leaf_positions1 = get_leaf_positions(tree1)
leaf_positions2 = get_leaf_positions(tree2)

# Create a figure
fig, axes = plt.subplots(1, 2, figsize=(12, 8))

# Draw the first tree on the left
Phylo.draw(tree1, do_show=False, axes=axes[0])
axes[0].set_title("Tree 1")

# Draw the second tree on the right
Phylo.draw(tree2, do_show=False, axes=axes[1])
axes[1].set_title("Tree 2")

# Add connecting lines between matching leaf nodes
fig.canvas.draw()  # Ensure the figure is fully drawn before adding lines
for leaf_name in leaf_positions1.keys():
    if leaf_name in leaf_positions2:
        # Get the positions of the matching leaves
        x1, y1 = axes[0].transData.transform((0, leaf_positions1[leaf_name]))
        x2, y2 = axes[1].transData.transform((1, leaf_positions2[leaf_name]))

        # Convert back to figure coordinates
        x1_fig, y1_fig = fig.transFigure.inverted().transform((x1, y1))
        x2_fig, y2_fig = fig.transFigure.inverted().transform((x2, y2))

        # Draw a line connecting the matching leaves
        line = plt.Line2D([x1_fig, x2_fig], [y1_fig, y2_fig], transform=fig.transFigure, color="gray", linestyle="--", linewidth=0.5)
        fig.lines.append(line)

plt.tight_layout()
plt.show()