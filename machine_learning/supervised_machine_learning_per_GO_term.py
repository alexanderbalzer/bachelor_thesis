import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
from goatools.obo_parser import GODag
from goatools.base import download_go_basic_obo
import statsmodels.api as sm
from sklearn.preprocessing import StandardScaler
import re
from statsmodels.stats.outliers_influence import variance_inflation_factor
from adjustText import adjust_text
from statsmodels.stats.multitest import multipletests
from datetime import datetime
from tqdm import tqdm
from matplotlib.patches import Patch
import warnings
def tqdm_warning(message, category, filename, lineno, file=None, line=None):
    tqdm.write(warnings.formatwarning(message, category, filename, lineno, line))

warnings.showwarning = tqdm_warning
import sys

def tqdm_excepthook(exctype, value, traceback):
    tqdm.write(f"Uncaught exception: {value}")

sys.excepthook = tqdm_excepthook

def load_obo():
    '''
    load the latest obo file to exchange the GO terms with their aspects
    '''
    obo_path = "pipeline/go.obo"
    if not os.path.exists(obo_path):
        obo_path = download_go_basic_obo()  
    # Load the GO DAG
    go_dag = GODag(obo_path)
    return go_dag

def get_go_aspect(go_id, go_dag):
    """
    Get the aspect of a GO term.
    """
    if go_id in go_dag:
        go_term = go_dag[go_id]
        return go_term.name
    else:
        return go_id
    
def sanitize_filename(name):
    # Replace any character that is not alphanumeric, underscore, or slash with underscore
    return re.sub(r'[^A-Za-z0-9_\/]', '_', name)

def compute_vif(X):
    """Compute Variance Inflation Factor (VIF) for each feature in the DataFrame X."""
    vif_data = pd.DataFrame()
    vif_data["feature"] = X.columns
    vif_data["VIF"] = [variance_inflation_factor(X.values, i) for i in range(X.shape[1])]
    return vif_data

def drop_one_hot_collinear(X):
    """
    Drop one column from each set of one-hot encoded columns to avoid collinearity.
    Assumes columns are named like 'category_value'.
    """
    cols_to_drop = []
    prefix_groups = {}
    for col in X.columns:
        prefix = col.split('_')[0]
        prefix_groups.setdefault(prefix, []).append(col)
    for group in prefix_groups.values():
        if len(group) > 1:
            # Drop the first column in the group
            cols_to_drop.append(group[0])
    return X.drop(columns=cols_to_drop)

def run(name, go_dag, dir):

    # Set the working directory
    working_dir = os.path.dirname(dir + name + "/")

    # Read the feature matrix from the CSV file
    feature_matrix_path = working_dir + "/feature_matrix_with_go_terms.csv"
    df = pd.read_csv(feature_matrix_path, index_col=0)

    # Drop the "Sequence" column if it exists
    if "Sequence" in df.columns:
        df = df.drop(columns=["Sequence"])
    if "protein_id_human" in df.columns:
        df = df.drop(columns=["protein_id_human"])

    # Create output directory for results
    general_output_dir = os.path.join(working_dir, "go_term_rf_results_mito_and_ER_2")
    os.makedirs(general_output_dir, exist_ok=True)
    '''
    go_ids = [
            #Go terms for: ER, Golgi, Ribosome, Mitochondria, Nucleus, Lysosome, Cell membrane, Cytoplasm
            "GO:0005783",
            "GO:0005794",
            "GO:0005840",
            "GO:0005739_no_cleavable_mts",
            "GO:0005739_cleavable_mts",
            "GO:0005634",
            "GO:0005764",
            "GO:0005886",
            "GO:0005737",
            "Multiple",
            "no_GO_term_found"
        ]
    valid_go_terms = [
            "GO:0005740",
            "GO:0005743",
            "GO:0005741",
            "GO:0005758",
            "GO:0005759",
            "GO:0005746",
            "GO:0031966",
            "Multiple"
        ]
    #valid_go_terms = go_ids
    '''
    valid_go_terms = ['GO:0005783', 'GO:0005739', 'Multiple']
    for go_term in tqdm(valid_go_terms, position=2, desc=f"Processing GO terms for {name}", leave=False):
        output_dir_go = os.path.join(general_output_dir, sanitize_filename(get_go_aspect(go_term, go_dag)))
        os.makedirs(output_dir_go, exist_ok=True)
        # Filter for current GO term (binary classification: this term vs. all others)
        df_filtered = df.copy()
        if go_term == 'GO:0005739' and name == 'Homo_sapiens':
            nat_types = ['natB', 'natA', 'natC', 'natX', 'all']
        else:
            nat_types = ['all']
        
        for nat_type in nat_types:
            df_filtered = df.copy()
            if go_term == "GO:0005739" and nat_type != "all":
                second_aa_columns = [col for col in df_filtered.columns if col.startswith("Second_AA_")]
                df_filtered["Second_AA"] = df_filtered[second_aa_columns].idxmax(axis=1).str.replace("Second_AA_", "")
                df_filtered = df_filtered.drop(columns=second_aa_columns)
                natA = {'A', 'S', 'T', 'V', 'C', 'G'}
                natC = {'L', 'I', 'F', 'W'}
                natB = {'D', 'E', 'N', 'Q'}
                # Add a column for nat type based on Second_AA
                def get_nat_type(aa):
                    if aa in natA:
                        return 'natA'
                    elif aa in natC:
                        return 'natC'
                    elif aa in natB:
                        return 'natB'
                    else:
                        return 'natX'  # Unknown or other types
                df_filtered['Nat_Type'] = df_filtered['Second_AA'].apply(get_nat_type)
                # delete the second amino acid column
                df_filtered = df_filtered.drop(columns=["Second_AA"])
                copy_of_df_filtered = df_filtered.copy()
                df_filtered = copy_of_df_filtered.copy()
                output_dir = os.path.join(output_dir_go, nat_type)
                os.makedirs(output_dir, exist_ok=True)
                # remove the rows that dont contain the go_term or the go term cyto_nuclear
                df_filtered = df_filtered[(df_filtered['GO_Term'] == go_term) | (df_filtered["GO_Term"] == "cyto_nuclear")]
                df_filtered = df_filtered[(df_filtered['Nat_Type'] == nat_type) | (df_filtered["GO_Term"] == "cyto_nuclear")]
                df_filtered["GO_Term_Binary"] = (df_filtered["GO_Term"] == go_term).astype(int)
                df_filtered = df_filtered.drop(columns=["Nat_Type"])
            else:
                output_dir = os.path.join(output_dir_go)
                os.makedirs(output_dir, exist_ok=True)
                # remove the rows that dont contain the go_term
                df_filtered = df_filtered[(df_filtered['GO_Term'] == go_term) | (df_filtered["GO_Term"] == "cyto_nuclear")]
                df_filtered["GO_Term_Binary"] = (df_filtered["GO_Term"] == go_term).astype(int)
                df_filtered = df_filtered.drop(columns=["Second_AA_V"])
        
            X = df_filtered.drop(['Molecular Weight', "GO_Term", "GO_Term_Binary", "Leucine_and_Alanine_percentage", "Arginine_percentage", "Discrimination Factor"], axis=1)
            y = df_filtered["GO_Term_Binary"]

            # Count how many proteins have the current GO term
            proteins_with_go_term = df_filtered[(df_filtered['GO_Term'] == go_term)].shape[0]

            dropped_constant_columns = X.columns[X.nunique() <= 1].tolist()
            dropped_duplicate_columns = X.T[X.T.duplicated()].index.tolist()

            if dropped_constant_columns:
                print(f"Dropping constant columns for {go_term}: {dropped_constant_columns}")
            if dropped_duplicate_columns:
                print(f"Dropping duplicate columns for {go_term}: {dropped_duplicate_columns}")

            X = X.loc[:, X.nunique() > 1]  # Remove constant columns
            X = X.T.drop_duplicates().T    # Remove duplicate columns

            # Add intercept for statsmodels
            '''X = drop_one_hot_collinear(X)'''
            X = sm.add_constant(X)
            
            # Compute and save VIF before standardization (exclude intercept if present)
            vif_df = compute_vif(X.drop(columns=["const"], errors="ignore"))
            vif_df.to_csv(os.path.join(output_dir, f"{go_term}_vif.csv"), index=False)

            # Z-transformation (standardization) of features (excluding the intercept)
            scaler = StandardScaler()
            X_scaled = scaler.fit_transform(X.drop(columns=["const"], errors="ignore"))
            X_scaled = pd.DataFrame(X_scaled, columns=[col for col in X.columns if col != "const"], index=X.index)
            X_scaled = sm.add_constant(X_scaled)  # Add intercept back

            # Fit logistic regression
            model = sm.Logit(y, X_scaled)
            try:
                result = model.fit(disp=0)
            except Exception as e:
                print(f"Skipping {go_term} due to fitting error: {e}")
                continue
            #go_term = sanitize_filename(get_go_aspect(go_term, go_dag))
            # Get coefficients (excluding intercept)
            coefs = result.params.drop("const")
            feature_importance_df = pd.DataFrame({
                'Feature': coefs.index,
                'Coefficient': coefs.values,
                'AbsCoefficient': np.abs(coefs.values),
                'Significance': result.pvalues[coefs.index]
            }).sort_values(by='AbsCoefficient', ascending=False)

            # FDR correction (Benjamini-Hochberg)
            rejected, pvals_corrected, _, _ = multipletests(feature_importance_df['Significance'], alpha=0.05, method='fdr_bh')
            feature_importance_df['FDR'] = pvals_corrected
            feature_importance_df['FDR_significant'] = rejected
            feature_importance_df['n'] = proteins_with_go_term

            # Save features with their significance and FDR to file
            feature_importance_df.to_csv(os.path.join(output_dir, f"{go_term}_logreg_coefficients_{name}.csv"), index=False)
            logreg_coefficients_path = os.path.join('pipeline/output/output_20250603_145910_ml_all_organisms', f"logreg_coefficients")
            os.makedirs(logreg_coefficients_path, exist_ok=True)
            feature_importance_df.to_csv(os.path.join(logreg_coefficients_path, f"{go_term}_logreg_coefficients_{name}.csv"), index=False)
            logreg_coeff = feature_importance_df[['Feature', 'Coefficient', 'FDR_significant']]

            # Save the model summary to a text file
            with open(os.path.join(output_dir, f"{go_term}_logreg_summary.txt"), "w") as f:
                f.write(result.summary().as_text())

            # Save mle_retvals
            with open(os.path.join(output_dir, f"{go_term}_mle_retvals.txt"), "w") as f:
                for key, val in result.mle_retvals.items():
                    f.write(f"{key}: {val}\n")

            # Plot and save feature importances
            '''plt.figure(figsize=(10, 6))
            sns.barplot(x='Coefficient', y='Feature', data=feature_importance_df)
            plt.title(f"Logistic Regression Coefficients for GO Term: {go_term}")
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, f"{go_term}_logreg_coefficients.png"))
            plt.close()'''
            feature_importance_df = feature_importance_df.replace([np.inf, -np.inf], np.nan)
            feature_importance_df = feature_importance_df.dropna(subset=["Significance", "Coefficient"])
            feature_importance_df = feature_importance_df[feature_importance_df["Significance"] > 0]

            plt.figure(figsize=(10, 6))
            sns.barplot(x='Coefficient', y='Feature', data=feature_importance_df.head(10))
            plt.title(f"Top 10 Logistic Regression Coefficients for GO Term: {go_term}")
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, f"{go_term}_logreg_coefficients_top10.png"))
            plt.close()

            # plot a volcano plot
            plt.figure(figsize=(10, 6))
            sns.scatterplot(x='Coefficient', y=-np.log10(feature_importance_df['Significance']), data=feature_importance_df, color='grey')
            if go_term == "GO:0005739":
                plt.title(f"Volcano Plot for {name}, {go_term}, {get_go_aspect(go_term, go_dag)}, {nat_type}, n = {proteins_with_go_term}")
            else:
                plt.title(f"Volcano Plot for {name}, {go_term}, {get_go_aspect(go_term, go_dag)}, n = {proteins_with_go_term}")
            plt.xlabel("Coefficient")
            plt.ylabel("-log10(Significance)")
            plt.axhline(y=-np.log10(0.05), color='r', linestyle='--')
            plt.axvline(x=0, color='g', linestyle='--')
            # Highlight significant features
            significant_features = feature_importance_df[(feature_importance_df['FDR'] < 0.1)]
            plt.scatter(significant_features['Coefficient'], -np.log10(significant_features['Significance']), color='red', label='Significant')

            # Use adjustText for non-overlapping labels
            texts = []
            for _, row in significant_features.iterrows():
                texts.append(
                    plt.text(row['Coefficient'], -np.log10(row['Significance']), row['Feature'], fontsize=8, color='blue')
                )
            adjust_text(texts, arrowprops=dict(arrowstyle='->', color='black', lw=0.5))

            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, f"{go_term}_volcano_plot.pdf"))
            plt.close()

        if go_term == "GO:0005739":
            mito_df = logreg_coeff
        elif go_term == "GO:0005783":
            er_df = logreg_coeff
        elif go_term == "Multiple":
            both_df = logreg_coeff

    '''# === 2. Gruppenzuweisung ===
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
    plt.savefig(os.path.join(general_output_dir, f"compared_coefficients_mito_ER_{name}.pdf"))'''
        

    for root, dirs, files in os.walk(general_output_dir, topdown=False):
        for name in dirs:
            dir_path = os.path.join(root, name)
            if not os.listdir(dir_path):
                os.rmdir(dir_path)
                print(f"Deleted empty directory: {dir_path}")

if __name__ == "__main__":
    organism_names = [
    "Homo_sapiens","Mus_musculus", "Rattus_norvegicus", "Dario_rerio",
    "Caenorhabditis_elegans", "Drosophila_Melanogaster", "Arabidopsis_thaliana", 
    "Saccharomyces_cerevisiae"]
    working_dir = 'pipeline/output/output_20250616_165204'
    #  "Homo_sapiens_isoforms", 
    #organism_names = ["Homo_sapiens"]
    # Load the GO DAG
    go_dag = load_obo()
    for name in tqdm(organism_names, position=1, desc= "pipeline progress", leave=False):
        run(name, go_dag, dir=working_dir)

