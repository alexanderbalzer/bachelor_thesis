import numpy as np
from scipy.stats import pearsonr
import numpy as np
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import pdist, squareform
from scipy.spatial.distance import euclidean, cosine
from io import StringIO
from scipy.cluster.hierarchy import to_tree
import os

def calculate_similarity_matrix_with_pearson(data):
    # Load the numpy array from the specified file
    n = len(data)
    similarity_matrix = np.zeros((n, n))

    for i in range(n):
        for j in range(n):
            if i == j:
                similarity_matrix[i][j] = 1.0  # Correlation with itself is 1
            elif i < j:  # Compute only for upper triangle
                corr, _ = pearsonr(data[i], data[j])
                similarity_matrix[i][j] = corr
                similarity_matrix[j][i] = corr  # Symmetric matrix

    return similarity_matrix

def calculate_similarity_matrix_with_euclidean(data):
    n = len(data)
    similarity_matrix = np.zeros((n, n))

    for i in range(n):
        for j in range(n):
            if i == j:
                similarity_matrix[i][j] = 1.0  # Similarity with itself is 1
            elif i < j:  # Compute only for upper triangle
                dist = euclidean(data[i], data[j])
                similarity = 1 / (1 + dist)  # Convert distance to similarity
                similarity_matrix[i][j] = similarity
                similarity_matrix[j][i] = similarity  # Symmetric matrix

    return similarity_matrix

def calculate_similarity_matrix_with_cosine(data):
    n = len(data)
    similarity_matrix = np.zeros((n, n))

    for i in range(n):
        for j in range(n):
            if i == j:
                similarity_matrix[i][j] = 1.0  # Similarity with itself is 1
            elif i < j:  # Compute only for upper triangle
                similarity = 1 - cosine(data[i], data[j])  # Cosine similarity
                similarity_matrix[i][j] = similarity
                similarity_matrix[j][i] = similarity  # Symmetric matrix

    return similarity_matrix

def convert_to_newick(linkage_matrix, labels):
    def build_newick(node, parent_dist, leaf_names):
        if node.is_leaf():
            return f"{leaf_names[node.id]}:{parent_dist - node.dist:.6f}"
        else:
            left = build_newick(node.left, node.dist, leaf_names)
            right = build_newick(node.right, node.dist, leaf_names)
        return f"({left},{right}):{parent_dist - node.dist:.6f}"

    tree = to_tree(linkage_matrix, rd=False)
    return build_newick(tree, tree.dist, labels) + ";"


def run(organism_names, cache_dir, output_dir, phylo_tree_method, phylo_tree_algorithm, save_newick):
    method = phylo_tree_method
    algorythm = phylo_tree_algorithm
    # Read the organisms list from the file
    labels = organism_names
    data = np.load(os.path.join(cache_dir, "phyl_tree_array.npy"))
    data = np.array(data)


    if method == "pearson":
        similarity_matrix = calculate_similarity_matrix_with_pearson(data)
    elif method == "euclidean":
        similarity_matrix = calculate_similarity_matrix_with_euclidean(data)
    elif method == "cosine":
        similarity_matrix = calculate_similarity_matrix_with_cosine(data)



    # Convert the similarity matrix to a distance matrix
    distance_matrix = 1 - similarity_matrix

    # Perform hierarchical clustering using UPGMA (Unweighted Pair Group Method with Arithmetic Mean)
    linkage_matrix = linkage(squareform(distance_matrix), method='average')


    newick_string = convert_to_newick(linkage_matrix, labels)


    # Save the Newick string to a file
    if save_newick:
        with open(os.path.join(output_dir, "phylogenetic_tree.newick"), "w") as file:
            file.write(newick_string)
    

    # Plot the phylogenetic tree
    plt.figure(figsize=(10, 7))
    dendrogram(linkage_matrix, labels= labels, orientation = 'left', leaf_font_size=10)
    plt.title("Phylogenetic Tree (UPGMA) - Pearson Correlation")
    plt.xlabel("Samples")
    plt.ylabel("Distance")
    plt.show()
    # Save the plot to a file
    plt.savefig(os.path.join(output_dir, "phylogenetic_tree.png"), dpi=300)

    # delete the cache
    os.remove(os.path.join(cache_dir, "phyl_tree_array.npy"))
    return

if __name__ == "__main__":
    organism_names = ["Organism1", "Organism2", "Organism3"]  # Example organism names
    cache_dir = "Phylogeny/cache"  # Example cache directory
    output_dir = "Phylogeny/output"  # Example output directory
    phylo_tree_method = "pearson"  # Example method
    phylo_tree_algorithm = "UPGMA"  # Example algorithm
    save_newick = True  # Example flag to save Newick format

    run(organism_names, cache_dir, output_dir, phylo_tree_method, phylo_tree_algorithm, save_newick)