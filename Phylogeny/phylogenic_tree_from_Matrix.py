import numpy as np
from scipy.stats import pearsonr
import numpy as np
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import pdist, squareform
from scipy.spatial.distance import euclidean, cosine

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


method = "pearson"  # Choose the method: "pearson", "euclidean", or "cosine"
# Read the organisms list from the file
with open("organisms_list.txt", "r") as file:
    labels = [line.strip() for line in file.readlines()]

data = np.load("subset_array.npy")
data = np.array(data)


if method == "pearson":
    similarity_matrix = calculate_similarity_matrix_with_pearson(data)
elif method == "euclidean":
    similarity_matrix = calculate_similarity_matrix_with_euclidean(data)
elif method == "cosine":
    similarity_matrix = calculate_similarity_matrix_with_cosine(data)
print("Similarity Matrix:")
print(similarity_matrix)



# Convert the similarity matrix to a distance matrix
distance_matrix = 1 - similarity_matrix

# Perform hierarchical clustering using UPGMA (Unweighted Pair Group Method with Arithmetic Mean)
linkage_matrix = linkage(squareform(distance_matrix), method='average')

# Plot the phylogenetic tree
plt.figure(figsize=(10, 7))
dendrogram(linkage_matrix, labels= labels, orientation = 'left', leaf_font_size=10)
plt.title("Phylogenetic Tree (UPGMA) - Pearson Correlation")
plt.xlabel("Samples")
plt.ylabel("Distance")
plt.show()
