import dendropy
from dendropy.interop import raxml
from matplotlib import gridspec

import matplotlib.pyplot as plt

def load_tree(file_path):
    """
    Load a phylogenetic tree from a Newick file.
    """
    return dendropy.Tree.get(path=file_path, schema="newick")

def plot_tanglegram(tree1, tree2):
    """
    Plot a tanglegram to compare two phylogenetic trees.
    """
    fig = plt.figure(figsize=(10, 6))
    gs = gridspec.GridSpec(1, 2, width_ratios=[1, 1])

    # Plot first tree
    ax1 = plt.subplot(gs[0])
    tree1.plot(ax=ax1, show_leaf_labels=True)
    ax1.set_title("Tree 1")

    # Plot second tree
    ax2 = plt.subplot(gs[1])
    tree2.plot(ax=ax2, show_leaf_labels=True)
    ax2.set_title("Tree 2")

    # Draw connecting lines
    for leaf1, leaf2 in zip(tree1.leaf_nodes(), tree2.leaf_nodes()):
        x1, y1 = ax1.transData.transform((0, leaf1.y))
        x2, y2 = ax2.transData.transform((1, leaf2.y))
        plt.plot([x1, x2], [y1, y2], color="gray", linestyle="--", linewidth=0.5)

    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    # Paths to the Newick tree files
    tree1_path = "path/to/tree1.newick"  # Replace with the output of phylogenic_tree_from_matrix.py
    tree2_path = "path/to/tree2.newick"  # Replace with the output of phylogenic_tree_from_aligned_fasta.py

    # Load the trees
    tree1 = load_tree(tree1_path)
    tree2 = load_tree(tree2_path)

    # Plot the tanglegram
    plot_tanglegram(tree1, tree2)