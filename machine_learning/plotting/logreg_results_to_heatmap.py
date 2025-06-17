import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os
from tqdm import tqdm

def format_species_name(name: str) -> str:
    # Teile den Namen anhand des Unterstrichs
    parts = name.split("_")
    if name == "Homo_sapiens_isoforms":
        return "H. Sapiens with isoforms"
    if len(parts) != 2:
        raise ValueError("Name muss genau ein Unterstrich enthalten (Gattung_Art)")
    genus, species = parts
    # Kürze den Gattungsnamen auf den ersten Buchstaben + Punkt
    short_genus = genus[0] + "."
    
    # Setze alles in kursiv (z. B. für Markdown oder HTML)
    formatted = f"{short_genus} {species}"
    return formatted

def run(organism_names, go_term, dir):
    combined_df = pd.DataFrame()
    amount_of_proteins = {}
    for name in organism_names:
        # Load the feature importance data
        logreg_coefficients_path = os.path.join(dir, f"logreg_coefficients")
        feature_importance_df_path = os.path.join(logreg_coefficients_path, f"{go_term}_logreg_coefficients_{name}.csv")
        feature_importance_df = pd.read_csv(feature_importance_df_path)
        
        # Keep only the columns 'Coefficient' and 'FDR_significant'
        feature_importance_df = feature_importance_df[['Feature','Coefficient', 'FDR_significant', 'n']]
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
        amount_of_proteins[name] = feature_importance_df.loc['Isoelectric Point', 'n']
        if 'mitofates cleavage probability' in feature_importance_df.index:
            # Remove the 'mitofates cleavage probability' and 'MPP cleavage position' rows
            feature_importance_df = feature_importance_df.drop(index=['mitofates cleavage probability', 'MPP cleavage position'])
        elif 'MPP cleavage position' in feature_importance_df.index:
            feature_importance_df = feature_importance_df.drop(index=['MPP cleavage position'])

        # Combine the dataframes
        combined_df = pd.concat([combined_df, feature_importance_df], ignore_index=False)
    

    # Format organism names for better readability and add amount of proteins
    combined_df['Organism'] = combined_df['Organism'].apply(lambda x: f"{format_species_name(x)} (n = {amount_of_proteins[x]})")
    # Pivot: rows=Organism, columns=Feature, values=Coefficient
    pivot_df = combined_df.pivot_table(
        index='Feature', columns='Organism', values='Coefficient', aggfunc='mean'
    )

    # Create a mask for significant values (same shape as pivot_df)
    significance_mask = combined_df.pivot_table(
        index='Feature', columns='Organism', values='FDR_significant', aggfunc=lambda x: any(x == True)
    )

    # Nur Features behalten, die mindestens einen signifikanten Wert haben
    rows_with_significance = significance_mask.any(axis=1)
    pivot_df = pivot_df.loc[rows_with_significance]
    significance_mask = significance_mask.loc[rows_with_significance]

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
        cbar_kws={'orientation': 'vertical', 'label': 'Coefficient', 'ticks': np.arange(-0.8, 0.9, 0.8)},
        linewidths=0.5,
        cbar_pos=(0.85, 0.83, 0.03, 0.15),  # [left, bottom, width, height]
        row_cluster=True,  # Cluster rows (features)
        col_cluster=True,  # Cluster columns (organisms),
    )

    # name the axes
    g.ax_heatmap.set_xlabel("")
    g.ax_heatmap.set_ylabel("")
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=45, ha='right')

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
    plot_folder = os.path.join(working_dir, "plots")
    os.makedirs(plot_folder, exist_ok=True)
    g.savefig(os.path.join(plot_folder, f"{go_term}_clustermap.pdf"), dpi=300)

if __name__ == "__main__":
    # Define the organism names and GO term
    organism_names = [
    "Homo_sapiens","Mus_musculus", "Rattus_norvegicus", "Danio_rerio",
    "Caenorhabditis_elegans", "Drosophila_Melanogaster", "Arabidopsis_thaliana", 
    "Saccharomyces_cerevisiae"]
    #go_term = 'GO:0005739' # one of GO:0005739 GO:0005783 Multiple
    go_terms = ['GO:0005739','GO:0005783', 'Multiple']
    working_dir = 'pipeline/output/output_20250617_151116_latest_ML'
    for go_term in go_terms:
        run(organism_names=organism_names, go_term=go_term, dir=working_dir)
