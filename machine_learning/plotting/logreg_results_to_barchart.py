import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os
from tqdm import tqdm
from matplotlib.patches import Patch


def run(organism_names, go_terms):
    for name in organism_names:
        for go_term in go_terms:
            # Load the feature importance data
            logreg_coefficients_path = os.path.join('pipeline/output/output_20250603_145910_ml_all_organisms', f"logreg_coefficients")
            feature_importance_df_path = os.path.join(logreg_coefficients_path, f"{go_term}_logreg_coefficients_{name}.csv")
            feature_importance_df = pd.read_csv(feature_importance_df_path)
            # Keep only the columns 'Coefficient' and 'FDR_significant'
            feature_importance_df = feature_importance_df[['Feature', 'Coefficient', 'FDR_significant']]
            feature_importance_df['Feature'] = feature_importance_df['Feature'].str.replace('_', ' ')
            # Sort the index alphabetically
            if go_term == "GO:0005739":
                mito_df = feature_importance_df
            elif go_term == "GO:0005783":
                er_df = feature_importance_df
            elif go_term == "Multiple":
                both_df = feature_importance_df

            # === 2. Gruppenzuweisung ===
        mito_df['Group'] = 'Mitochondrion (GO:0005739)'
        er_df['Group'] = 'ER (GO:0005783)'
        both_df['Group'] = 'ER+Mito'

        # === 3. Kombinieren und nur gemeinsame Features behalten ===
        combined_df = pd.concat([mito_df, er_df, both_df], ignore_index=True)
        common_features = set(mito_df['Feature']) & set(er_df['Feature']) & set(both_df['Feature'])
        combined_df = combined_df[combined_df['Feature'].isin(common_features)]

        # === 4. Zwei DataFrames erstellen: signifikant und alle Werte ===
        full_pivot = combined_df.pivot(index='Feature', columns='Group', values='Coefficient')
        sig_df = combined_df.copy()
        sig_df.loc[~sig_df['FDR_significant'], 'Coefficient'] = np.nan
        sig_pivot = sig_df.pivot(index='Feature', columns='Group', values='Coefficient')

        # === 5. Plot erstellen ===
        fig, ax = plt.subplots(figsize=(8, 10))
        # Zuerst: alle (grau)
        full_pivot.plot(kind='barh', ax=ax, edgecolor='white', legend=False, hatch='/////', width=0.8)

        # Dann: signifikante Balken drüberlegen (farbig)
        sig_pivot.plot(kind='barh', ax=ax, edgecolor='white', legend=True, width=0.8)
        # Add a second legend for hatched and unhatched bars
        legend_elements = [
            Patch(facecolor='blue', edgecolor='black', label='ER (GO:0005783)'),
            Patch(facecolor='orange', edgecolor='black', label='ER+Mito'),
            Patch(facecolor='green', edgecolor='black', label='Mitochondrion (GO:0005739)'),
            Patch(facecolor='none', edgecolor='none', label=''),
            Patch(facecolor='gray', edgecolor='white', hatch='/////', label='Not Significant'),
            Patch(facecolor='gray', edgecolor='black', label='Significant')
        ]
        ax.legend(handles=legend_elements, loc='upper right', fontsize='small')

        # === 6. Plot verschönern ===
        plt.xlabel("Coefficient")
        plt.ylabel("Feature")
        plt.axvline(x=0, color='black', linestyle='-', linewidth=0.8)
        plt.xticks(rotation=0)
        plt.tight_layout()
        plt.savefig(os.path.join(logreg_coefficients_path, f"{go_term}_logreg_coefficients_{name}.pdf"), dpi=300)

if __name__ == "__main__":
    # Define the organism names and GO term
    organism_names = [
    "Homo_sapiens", "Mus_musculus", "Dario_rerio", "Daphnia_magna", 
    "Caenorhabditis_elegans", "Drosophila_Melanogaster", "Arabidopsis_thaliana", 
    "Physcomitrium_patens", "Chlamydomonas_reinhardtii", 
    "Candida_glabrata", "Saccharomyces_cerevisiae", "Zygosaccharomyces_rouxii"]
    organism_names = ['Homo_sapiens']
    go_terms = ['GO:0005739','GO:0005783', 'Multiple']
    run(organism_names=organism_names, go_terms=go_terms)
