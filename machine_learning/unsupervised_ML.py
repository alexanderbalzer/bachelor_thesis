import numpy as np
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button
import seaborn as sns
import pandas as pd
import os
from sklearn.cluster import KMeans

name = "Homo_sapiens"  # Example organism name
# Set the working directory
working_dir = os.path.dirname("pipeline/output/output_20250515_105213/" + name + "/")
feature_matrix_path = working_dir + "/feature_matrix.csv"
# Read the feature matrix from the CSV file
feature_matrix = pd.read_csv(feature_matrix_path, index_col=0)

# Save the original cleavable_mts info
cleavable_values = feature_matrix["cleavable_mts"].values
feature_matrix = feature_matrix.drop(columns=["cleavable_mts"])

# Drop non-numeric columns if any
feature_matrix_numeric = feature_matrix.select_dtypes(include=[np.number])

scaler = StandardScaler()
data_scaled = scaler.fit_transform(feature_matrix_numeric)

# Perform PCA (with 2 components)
pca = PCA(n_components=2)
principal_components = pca.fit_transform(data_scaled)
principal_df = pd.DataFrame(principal_components, columns=['PC1', 'PC2'])

# Explained variance (How much variance is explained by each principal component?)
explained_variance = pca.explained_variance_ratio_

# Print explained variance
print("Explained variance of components:\n", explained_variance)

# After fitting PCA
feature_names = feature_matrix_numeric.columns

# Identify which features contribute most to each principal component
feature_contributions = pd.DataFrame(
    pca.components_,
    columns=feature_names,
    index=['PC1', 'PC2']
)

# Print the contributions of features to each principal component
print("Feature contributions to principal components:")
print(feature_contributions)

# Find the top contributing features for each principal component
top_features_pc1 = feature_contributions.loc['PC1'].nlargest(5)
top_features_pc2 = feature_contributions.loc['PC2'].nlargest(5)

print("\nTop contributing features for PC1:")
print(top_features_pc1)

print("\nTop contributing features for PC2:")
print(top_features_pc2)
# --- Interactive plot with slider and save button ---
fig, ax = plt.subplots(figsize=(8, 6))
plt.subplots_adjust(bottom=0.3)  # Make space for the slider and button

# Initial threshold
init_threshold = np.median(cleavable_values)

# Function to update plot
def update(threshold):
    labels = np.where(cleavable_values > threshold, "yes", "no")
    ax.clear()
    sns.scatterplot(
        x=principal_df["PC1"], y=principal_df["PC2"],
        hue=labels, palette={"yes": "tab:blue", "no": "tab:orange"},
        edgecolor="black", ax=ax
    )
    ax.set_title('PCA: Visualization of the First 2 Principal Components')
    ax.set_xlabel('PC1 (First Principal Component)')
    ax.set_ylabel('PC2 (Second Principal Component)')
    ax.legend(title="Cleavable MTS")
    ax.grid(True)
    ax.text(
        0.5, 0.5,
        f'Explained Variance:\nPC1: {pca.explained_variance_ratio_[0]:.2f}\nPC2: {pca.explained_variance_ratio_[1]:.2f}',
        horizontalalignment='center', verticalalignment='center', fontsize=12,
        bbox=dict(facecolor='white', alpha=0.7),
        transform=ax.transAxes
    )
    fig.canvas.draw_idle()

# Initial plot
update(init_threshold)

# Slider
axcolor = 'lightgoldenrodyellow'
ax_thresh = plt.axes([0.2, 0.12, 0.6, 0.03], facecolor=axcolor)
slider = Slider(ax_thresh, 'Threshold', float(np.min(cleavable_values)), float(np.max(cleavable_values)), valinit=init_threshold)

slider.on_changed(update)

# Save button
ax_save = plt.axes([0.8, 0.02, 0.15, 0.06])
button = Button(ax_save, 'Save Plot', color='lightblue', hovercolor='skyblue')

def save_plot(event):
    threshold = slider.val
    filename = f"PCA_plot_threshold_{threshold:.2f}.pdf"
    file = os.path.join(working_dir, filename)
    plt.savefig(file)
    print(f"Plot saved as {filename}")

button.on_clicked(save_plot)

plt.show()


# Perform K-Means clustering on the first two principal components
kmeans = KMeans(n_clusters=2, random_state=42)
clusters = kmeans.fit_predict(principal_df)

# Add cluster labels to the DataFrame
principal_df['Cluster'] = clusters

# Print cluster centers
print("Cluster centers (in PC1 and PC2 space):")
print(kmeans.cluster_centers_)

# Visualize the clusters
fig, ax = plt.subplots(figsize=(8, 6))
sns.scatterplot(
    x=principal_df["PC1"], y=principal_df["PC2"],
    hue=principal_df["Cluster"], palette="viridis", edgecolor="black", ax=ax
)

# Plot the cluster centers
cluster_centers = kmeans.cluster_centers_
ax.scatter(
    cluster_centers[:, 0], cluster_centers[:, 1],
    s=200, c='red', marker='X', label='Cluster Centers', edgecolor='black'
)

ax.legend(title="Cluster")
ax.set_title('K-Means Clustering on PCA Components')
ax.set_xlabel('PC1 (First Principal Component)')
ax.set_ylabel('PC2 (Second Principal Component)')
ax.legend(title="Cluster")
ax.grid(True)
plt.show()