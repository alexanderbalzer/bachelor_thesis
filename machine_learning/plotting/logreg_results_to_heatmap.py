import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os
from tqdm import tqdm


def run(organism_names, go_term):
    combined_df = pd.DataFrame()
    for name in organism_names:
        # Load the feature importance data
        logreg_coefficients_path = os.path.join('pipeline/output/output_20250603_145910_ml_all_organisms', f"logreg_coefficients")
        feature_importance_df_path = os.path.join(logreg_coefficients_path, f"{go_term}_logreg_coefficients_{name}.csv")
        feature_importance_df = pd.read_csv(feature_importance_df_path)
        
        # Keep only the columns 'Coefficient' and 'FDR_significant'
        feature_importance_df = feature_importance_df[['Feature','Coefficient', 'FDR_significant']]
        feature_importance_df['Feature'] = feature_importance_df['Feature'].str.replace('_', ' ')
        # make the 'Feature' column the index
        feature_importance_df.set_index('Feature', inplace=True)
        # Sort the index alphabetically
        feature_importance_df.sort_index(inplace=True)
        # Add a column for the organism name
        feature_importance_df['Organism'] = name
        # Add a column for the coefficient
        feature_importance_df['Coefficient'] = feature_importance_df['Coefficient'].astype(float)
        # Add a column for significance
        feature_importance_df['FDR_significant'] = feature_importance_df['FDR_significant'].astype(str).str.lower() == 'true'
        # Combine the dataframes
        combined_df = pd.concat([combined_df, feature_importance_df], ignore_index=False)
    # Pivot: rows=Organism, columns=Feature, values=Coefficient
    pivot_df = combined_df.pivot_table(
        index='Feature', columns='Organism', values='Coefficient', aggfunc='mean'
    )

    # Create a mask for significant values (same shape as pivot_df)
    significance_mask = combined_df.pivot_table(
        index='Feature', columns='Organism', values='FDR_significant', aggfunc=lambda x: any(x == True)
    )

    '''plt.figure(figsize=(len(pivot_df.columns)*0.7+4, len(pivot_df.index)*0.5+4))
    ax = sns.heatmap(pivot_df, cmap='vlag', annot=False, cbar=True)

    # Add asterisks for significant values
    for y, organism in enumerate(pivot_df.index):
        for x, feature in enumerate(pivot_df.columns):
            if significance_mask.loc[organism, feature]:
                ax.text(x + 0.5, y + 0.5, '*', ha='center', va='center', color='black', fontsize=12)

    plt.xlabel("Feature")
    plt.ylabel("Organism")
    plt.title(f"LogReg Coefficient Heatmap for GO Term: {go_term}")
    plt.tight_layout()
    plt.savefig(os.path.join(logreg_coefficients_path, f"{go_term}_heatmap.png"), dpi=300)
    plt.show()'''

    # Create a clustermap (features=y, organisms=x, coefficients=values)
    g = sns.clustermap(
        pivot_df, 
        cmap='RdBu_r',
        vmin=-0.8, vmax=0.8,  # <-- set colorbar range
        figsize=(len(pivot_df.columns)*0.7+4, len(pivot_df.index)*0.5+4),
        cbar_kws={'label': 'Coefficient'},
        linewidths=0.5
    )
    # Add asterisks for significant values

    # Reorder index and columns based on clustermap dendrogram
    row_order = g.dendrogram_row.reordered_ind
    col_order = g.dendrogram_col.reordered_ind

    features_ordered = pivot_df.index[row_order]
    organisms_ordered = pivot_df.columns[col_order]

    # Add asterisks at the correct clustered positions
    for y, feature in enumerate(features_ordered):
        for x, organism in enumerate(organisms_ordered):
            if significance_mask.loc[feature, organism]:
                g.ax_heatmap.text(x + 0.5, y + 0.5, '*', ha='center', va='center', color='black', fontsize=12)


    g.ax_heatmap.set_xlabel("Organism")
    g.ax_heatmap.set_ylabel("Feature")

    g.savefig(os.path.join(logreg_coefficients_path, f"{go_term}_clustermap.pdf"), dpi=300)

if __name__ == "__main__":
    # Define the organism names and GO term
    organism_names = [
    "Homo_sapiens", "Mus_musculus", "Dario_rerio", "Daphnia_magna", 
    "Caenorhabditis_elegans", "Drosophila_Melanogaster", "Arabidopsis_thaliana", 
    "Physcomitrium_patens", "Chlamydomonas_reinhardtii", 
    "Candida_glabrata", "Saccharomyces_cerevisiae", "Zygosaccharomyces_rouxii"]
    #go_term = 'GO:0005739' # one of GO:0005739 GO:0005783 Multiple
    go_terms = ['GO:0005739','GO:0005783', 'Multiple']
    for go_term in go_terms:
        run(organism_names=organism_names, go_term=go_term)
