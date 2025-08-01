import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os
from tqdm import tqdm
from matplotlib.patches import Patch


def run(organism_names, go_terms, working_dir):
    for name in organism_names:
        for go_term in go_terms:
            # Load the feature importance data
            logreg_coefficients_path = os.path.join(working_dir, f"logreg_coefficients")
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
        #both_df['Group'] = 'ER+Mito'

        # === 3. Kombinieren und nur gemeinsame Features behalten ===
        combined_df = pd.concat([mito_df, er_df], ignore_index=True)
        common_features = set(mito_df['Feature']) & set(er_df['Feature'])
        combined_df = combined_df[combined_df['Feature'].isin(common_features)]
        # Set the first row as the header

        # === 4. Zwei DataFrames erstellen: signifikant und alle Werte ===
        full_pivot = combined_df.pivot(index='Feature', columns='Group', values='Coefficient')
        sig_df = combined_df.copy()
        sig_df.loc[~sig_df['FDR_significant'], 'Coefficient'] = np.nan
        sig_pivot = sig_df.pivot(index='Feature', columns='Group', values='Coefficient')

        # === 5. Plot erstellen ===
        color_map = {
            'ER (GO:0005783)': 'blue',
            'Mitochondrion (GO:0005739)': 'green'
        }
        fig, ax = plt.subplots(figsize=(12, 6))
        # First: all (gray)
        full_pivot.plot(kind='bar', ax=ax, edgecolor='white', legend=False, hatch='/////', width=0.8, color=[color_map.get(col, 'gray') for col in full_pivot.columns])

        # Then: overlay significant bars (colored)
        sig_pivot.plot(kind='bar', ax=ax, edgecolor='white', legend=True, width=0.8, color=[color_map.get(col, 'gray') for col in sig_pivot.columns])

        # Add a second legend for hatched and unhatched bars
        legend_elements = [
            Patch(facecolor='blue', edgecolor='black', label='ER (GO:0005783)'),
            Patch(facecolor='green', edgecolor='black', label='Mitochondrion (GO:0005739)'),
            Patch(facecolor='none', edgecolor='none', label=''),
            Patch(facecolor='gray', edgecolor='white', hatch='/////', label='Not Significant'),
            Patch(facecolor='gray', edgecolor='black', label='Significant')
        ]
        ax.legend(handles=legend_elements, loc='upper right', fontsize='small')

        # === 6. Plot verschönern ===
        plt.xlabel("Feature")
        plt.ylabel("Coefficient")
        plt.axhline(y=0, color='black', linestyle='-', linewidth=0.8)
        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()
        barchart_folder = os.path.join(working_dir, "plots")
        os.makedirs(barchart_folder, exist_ok=True)
        plt.savefig(os.path.join(barchart_folder, f"{go_term}_logreg_coefficients_{name}.pdf"), dpi=300)

if __name__ == "__main__":
    # Get the absolute path to the folder this script is in
    script_dir = os.path.dirname(os.path.abspath(__file__))

    # Build the relative path to the data file
    pipeline_output_path = os.path.join(script_dir, '..', 'pipeline', 'output')
    pipeline_input_path = os.path.join(script_dir, '..', 'pipeline', 'input')

    dirs = [d for d in os.listdir(pipeline_output_path) if os.path.isdir(os.path.join(pipeline_output_path, d))]

    # Sort alphabetically and get the last one
    if dirs:
        last_dir = sorted(dirs)[-1]
        last_output_dir = os.path.join(pipeline_output_path, last_dir)
        print("Last directory:", last_output_dir)
    else:
        print("No directories found.")

    # Normalize the path (optional but good practice)
    working_dir = os.path.normpath(last_output_dir)
    input_dir = pipeline_input_path
    # Read organism names from FASTA file names
    fasta_files = [f for f in os.listdir(input_dir) if f.endswith(".fasta")]
    organism_names = [os.path.splitext(f)[0] for f in fasta_files]
    go_terms = ['GO:0005739','GO:0005783']
    run(organism_names=organism_names, go_terms=go_terms, working_dir=working_dir)
