import numpy as np
from Bio import SeqIO
from scipy.stats import pearsonr
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import to_tree


def calculate_similarity_matrix(fasta_file):
    # Parse sequences from the FASTA file
    sequences = []
    ids = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        ids.append(" ".join(record.description.split()[1:3]))
        sequences.append(str(record.seq))
    
    # Ensure all sequences are of the same length
    seq_length = len(sequences[0])
    if not all(len(seq) == seq_length for seq in sequences):
        raise ValueError("All sequences must be of the same length.")
    
    # Convert sequences to numerical representation
    def seq_to_numeric(seq):
        mapping = {'a': 1, 'c': 2, 'g': 3, 't': 4, '-': 0}
        return [mapping.get(base, 0) for base in seq]  # Default to 0 for unknown bases
    
    numeric_sequences = [seq_to_numeric(seq) for seq in sequences]
    print(numeric_sequences)
    
    # Calculate the Pearson similarity matrix
    num_seqs = len(numeric_sequences)
    similarity_matrix = np.zeros((num_seqs, num_seqs))
    
    for i in range(num_seqs):
        for j in range(num_seqs):
            if i <= j:  # Compute only for upper triangle and diagonal
                corr, _ = pearsonr(numeric_sequences[i], numeric_sequences[j])
                similarity_matrix[i, j] = corr
                similarity_matrix[j, i] = corr  # Symmetric matrix
    
    return ids, similarity_matrix

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



def main():
    fasta_file = "aligned_sequences.fasta"  # Replace with your FASTA file path
    ids, similarity_matrix = calculate_similarity_matrix(fasta_file)
    


    # Print the similarity matrix
    print("Similarity Matrix:")
    print("\t" + "\t".join(ids))
    for i, row in enumerate(similarity_matrix):
        print(ids[i] + "\t" + "\t".join(f"{value:.3f}" for value in row))
    # Convert the similarity matrix to a distance matrix
    import matplotlib.pyplot as plt


    # Convert the similarity matrix to a distance matrix
    distance_matrix = 1 - similarity_matrix
    # Convert the square distance matrix to a condensed form
    condensed_distance_matrix = squareform(distance_matrix, checks=False)
    linkage_matrix = linkage(condensed_distance_matrix, method='average')

    newick_string = convert_to_newick(linkage_matrix, ids)


    # Save the Newick string to a file
    with open("Phylogeny/output/phylogenetic_tree_from_16s.newick", "w") as file:
        file.write(newick_string)

    # Plot the phylogenetic tree
    plt.figure(figsize=(10, 7))
    dendrogram(linkage_matrix, labels= ids, orientation = 'left', leaf_font_size=10)
    plt.title("Phylogenetic Tree (UPGMA) - Pearson Correlation")
    plt.xlabel("Samples")
    plt.ylabel("Distance")
    plt.show()


if __name__ == "__main__":
    main()

