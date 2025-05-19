import numpy as np
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import os
import umap

name = "Homo_sapiens"  # Example organism name
# Set the working directory
#working_dir = os.path.dirname("pipeline/output/output_20250514_134354/" + name + "/")
working_dir = os.path.dirname("pipeline/output/output_20250515_105213/" + name + "/")

feature_matrix_path = working_dir + "/feature_matrix_with_go_terms.csv"
# Read the feature matrix from the CSV file
feature_matrix = pd.read_csv(feature_matrix_path, index_col=0)
print(feature_matrix.head(10))
# Save the original cleavable_mts info (optional, not needed for GO_term coloring)
# cleavable_values = feature_matrix["cleavable_mts"].values
#feature_matrix = feature_matrix.drop(columns=["cleavable_mts"])
# Remove rows with duplicate protein IDs
feature_matrix = feature_matrix[~feature_matrix.index.duplicated(keep=False)]

# Save the GO_term column for later use
go_term_values = feature_matrix["GO_Term"].values
# Drop non-numeric columns if any, but keep GO_term for coloring
feature_matrix_numeric = feature_matrix.select_dtypes(include=[np.number])

scaler = StandardScaler()
data_scaled = scaler.fit_transform(feature_matrix_numeric)

# --- UMAP dimensionality reduction ---
umap_model = umap.UMAP(n_components=2, random_state=42)
umap_components = umap_model.fit_transform(data_scaled)
umap_df = pd.DataFrame(umap_components, columns=['UMAP1', 'UMAP2'])

# Add GO_term column for coloring
umap_df["GO_term"] = go_term_values

# --- Scatter plot colored by GO_term ---
plt.figure(figsize=(10, 7))
sns.scatterplot(
    x=umap_df["UMAP1"], y=umap_df["UMAP2"],
    hue=umap_df["GO_term"],
    palette="tab20",  # or another palette, depending on number of unique GO terms
    edgecolor="black"
)
plt.title('UMAP: Visualization colored by GO_term')
plt.xlabel('UMAP1')
plt.ylabel('UMAP2')
plt.legend(title="GO_term", bbox_to_anchor=(1.05, 1), loc='upper left')
plt.grid(True)
plt.tight_layout()
plt.show()