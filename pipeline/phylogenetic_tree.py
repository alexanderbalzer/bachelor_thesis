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
from utils import transform_labels_to_names
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, DistanceMatrix
from Bio import Phylo

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

    if algorythm == "UPGMA":
        # Perform hierarchical clustering using UPGMA (Unweighted Pair Group Method with Arithmetic Mean)
        linkage_matrix = linkage(squareform(distance_matrix), method='average')


        newick_string = convert_to_newick(linkage_matrix, labels)


        # Save the Newick string to a file
        if save_newick:
            with open(os.path.join(output_dir, "phylogenetic_tree.newick"), "w") as file:
                file.write(newick_string)


        # Plot the phylogenetic tree
        plt.figure(figsize=(10, 7))

        italic_labels = [f"$\\mathit{{{name.replace(' ', '\\ ')}}}$" for name in transform_labels_to_names(labels)]

        dendrogram(
            linkage_matrix, 
            labels= italic_labels, 
            orientation = 'left', 
            leaf_font_size=10
            )
        plt.title("Phylogenetic Tree (UPGMA) - Pearson Correlation")
        plt.xlabel("Samples")
        plt.ylabel("Distance")
        # Adjust the layout to prevent labels from being cut off
        plt.subplots_adjust(left=0.3)  # Increase the left margin
        plt.tight_layout()
        # Save the plot to a file
        plt.savefig(os.path.join(output_dir, "phylogenetic_tree.png"), dpi=300)

        # delete the cache
        os.remove(os.path.join(cache_dir, "phyl_tree_array.npy"))
        return
    
    elif algorythm == "nj":
        # Perform hierarchical clustering using Neighbor Joining (NJ)
        # Convert the similarity matrix to a distance matrix
        distance_matrix = [[1 - similarity_matrix[i][j] for j in range(len(similarity_matrix))] for i in range(len(similarity_matrix))]
        constructor = DistanceTreeConstructor()
        lower_triangle = []
        for i in range(len(distance_matrix)):
            lower_triangle.append(distance_matrix[i][:i + 1])  # Include only elements up to the diagonal
        italic_labels = [f"$\\mathit{{{name.replace(' ', '\\ ')}}}$" for name in transform_labels_to_names(labels)]
        distance_matrix = DistanceMatrix(names=labels, matrix=lower_triangle)
        tree = constructor.nj(distance_matrix)

        # Save the Newick string to a file
        if save_newick:
            with open(os.path.join(output_dir, "phylogenetic_tree.newick"), "w") as file:
                Phylo.write(tree, file, format="newick")

        distance_matrix = DistanceMatrix(names=italic_labels, matrix=lower_triangle)
        tree = constructor.nj(distance_matrix)
        # Plot the phylogenetic tree
        plt.figure(figsize=(10, 7))
        Phylo.draw(tree, do_show=False, label_func=lambda x: x.name)

        plt.title("Phylogenetic Tree")
        # Save the plot to a file
        plt.savefig(os.path.join(output_dir, "phylogenetic_tree.png"), dpi=300)

        # delete the cache
        os.remove(os.path.join(cache_dir, "phyl_tree_array.npy"))
        return


if __name__ == "__main__":
    # Example usage
    organism_names = ["Arabidopsis_thaliana",
"Caenorhabditis_elegans",
"Candida_glabrata",
"Clavispora_lusitaniae",
"Debaryomyces_hansenii",
"Drosophila_Melanogaster",
"Geotrichum_candidum",
"human",
"human_with_isoforms",
"Lachancea_thermotolerans",
"Mus_musculus",
"Physcomitrium_patens",
"Saccharomyces_cerevisiae",
"Scheffersomyces_stipitis",
"Schizosaccharomyces_pombe",
"Yarrowia_lipolytica",
"Zygosaccharomyces_rouxii"]
    cache_dir = "pipeline/cache/cache_20250414_224234/"
    output_dir = "pipeline/output/output_20250414_224234/"
    phylo_tree_method = "pearson"  # or "euclidean" or "cosine"
    phylo_tree_algorithm = "UPGMA"  # or "nj"
    save_newick = True

    run(organism_names, cache_dir, output_dir, phylo_tree_method, phylo_tree_algorithm, save_newick)