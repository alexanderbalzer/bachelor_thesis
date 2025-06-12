import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os
from tqdm import tqdm


def run(organism_names, go_term):
    for name in organism_names:
        # Load the feature importance data
        logreg_coefficients_path = os.path.join('pipeline/output/output_20250603_145910_ml_all_organisms', f"logreg_coefficients")
        feature_importance_df_path = os.path.join(logreg_coefficients_path, f"{go_term}_logreg_coefficients_{name}.csv")
        feature_importance_df = pd.read_csv(feature_importance_df_path)
        # Use the first row as column names and the first column as index
        feature_importance_df.columns = feature_importance_df.iloc[0]
        feature_importance_df = feature_importance_df[1:]
        feature_importance_df.set_index(feature_importance_df.columns[0], inplace=True)
        
        # Sort the index alphabetically
        feature_importance_df.sort_index(inplace=True)
        # Keep only the columns 'Coefficient' and 'FDR_significant'
        feature_importance_df = feature_importance_df[['Coefficient', 'FDR_significant']]

        # Append to a list to combine later
        if 'combined_df' not in locals():
            combined_df = feature_importance_df
        else:
            combined_df = pd.concat([combined_df, feature_importance_df])

        # Pivot the combined dataframe for clustermap
        pivot_df = combined_df.pivot_table(index=combined_df.index, values='Coefficient', aggfunc='mean')

        # Create a mask for significant values
        significance_mask = combined_df.pivot_table(index=combined_df.index, values='FDR_significant', aggfunc=lambda x: any(x == True))

        # Create a seaborn clustermap
        clustermap = sns.clustermap(pivot_df, cmap='vlag', metric='euclidean', method='average')

        # Add asterisks for significant values
        for (i, j), value in np.ndenumerate(significance_mask.values):
            if value:  # If True, add an asterisk
                clustermap.ax_heatmap.text(j + 0.5, i + 0.5, '*', ha='center', va='center', color='black')

        plt.title(f"Clustermap for GO Term: {go_term}")
        plt.show()

if __name__ == "__main__":
    # Define the organism names and GO term
    organism_names = [
    "Homo_sapiens", "Mus_musculus", "Dario_rerio", "Daphnia_magna", 
    "Caenorhabditis_elegans", "Drosophila_Melanogaster", "Arabidopsis_thaliana", 
    "Physcomitrium_patens", "Chlamydomonas_reinhardtii", 
    "Candida_glabrata", "Saccharomyces_cerevisiae", "Zygosaccharomyces_rouxii"]
    go_term = 'GO:0005739' # one of GO:0005739 GO:0005783 Multiple
    run(organism_names=organism_names, go_term=go_term)
