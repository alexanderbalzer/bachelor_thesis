import numpy as np
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import os
import umap
from sklearn.decomposition import PCA

name = "Homo_sapiens"  # Example organism name
# Set the working directory
#working_dir = os.path.dirname("pipeline/output/output_20250514_134354/" + name + "/")
working_dir = os.path.dirname("pipeline/output/output_20250519_142700_machine_learning_human/" + name + "/")

feature_matrix_path = working_dir + "/feature_matrix_with_go_terms.csv"
# Read the feature matrix from the CSV file
feature_matrix = pd.read_csv(feature_matrix_path, index_col=0)
# Save the original cleavable_mts info (optional, not needed for GO_term coloring)
# cleavable_values = feature_matrix["cleavable_mts"].values
#feature_matrix = feature_matrix.drop(columns=["cleavable_mts"])
# Remove rows with duplicate protein IDs
feature_matrix = feature_matrix[~feature_matrix.index.duplicated(keep=False)]

# Filter out rows with missing or empty GO_Term
feature_matrix = feature_matrix[feature_matrix["GO_Term"].notnull() & (feature_matrix["GO_Term"] != "")]

# replace GO:0005739_cleavable_mts with GO:0005739
feature_matrix["GO_Term"] = feature_matrix["GO_Term"].replace("GO:0005739_cleavable_mts", "GO:0005739")
feature_matrix["GO_Term"] = feature_matrix["GO_Term"].replace("GO:0005739_no_cleavable_mts", "GO:0005739")
feature_matrix["GO_Term"] = feature_matrix["GO_Term"].replace("no_GO_term_found", "")
feature_matrix["GO_Term"] = feature_matrix["GO_Term"].replace("Multiple", "")
feature_matrix["GO_Term"] = feature_matrix["GO_Term"].replace("GO:0005634", "")
feature_matrix["GO_Term"] = feature_matrix["GO_Term"].replace("GO:0005634", "")
feature_matrix["GO_Term"] = feature_matrix["GO_Term"].replace("GO:0005764", "GO:0005783")
feature_matrix["GO_Term"] = feature_matrix["GO_Term"].replace("GO:0005886", "GO:0005783")
feature_matrix["GO_Term"] = feature_matrix["GO_Term"].replace("GO:0005794", "GO:0005783")

# Filter out rows where GO_Term is empty
feature_matrix = feature_matrix[feature_matrix["GO_Term"] != ""]
# Save the GO_term column for later use
go_term_values = feature_matrix["GO_Term"].values

# Drop non-numeric columns if any, but keep GO_term for coloring
feature_matrix_numeric = feature_matrix.select_dtypes(include=[np.number])
# drop all columns except the specified ones
columns_to_keep = ["Hydrophobic Moment", "Hydrophobicity", "Isoelectric Point"]  # Replace with the actual column names you want to keep
feature_matrix_numeric = feature_matrix_numeric[columns_to_keep]

scaler = StandardScaler()
data_scaled = scaler.fit_transform(feature_matrix_numeric)

# --- PCA dimensionality reduction ---
pca = PCA(n_components=2, random_state=42)
pca_components = pca.fit_transform(data_scaled)
umap_df = pd.DataFrame(pca_components, columns=['PC1', 'PC2'])

# Add GO_term column for coloring
umap_df["GO_term"] = go_term_values

print(umap_df[['PC1', 'PC2', 'GO_term']])

# --- Scatter plot colored by GO_term ---
plt.figure(figsize=(10, 7))
sns.scatterplot(
    x=umap_df["PC1"], y=umap_df["PC2"],
    hue=umap_df["GO_term"],
    palette="tab20",
    edgecolor="black"
)
plt.title('PCA: Visualization colored by GO_term')
plt.xlabel('PC1')
plt.ylabel('PC2')
plt.legend(title="GO_term", bbox_to_anchor=(1.05, 1), loc='upper left')
# Save legend content to a text file
legend_labels = [text.get_text() for text in plt.gca().get_legend().get_texts()]
legend_file_path = os.path.join(working_dir, "legend_content.txt")
with open(legend_file_path, "w") as f:
    for label in legend_labels:
        f.write(label + "\n")
plt.grid(True)
plt.tight_layout(rect=[0, 0, 0.85, 1])
plt.show()