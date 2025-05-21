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

name = "Homo_sapiens"  # Example organism name


# Load the GO DAG
go_dag = load_obo()

# Set the working directory
working_dir = os.path.dirname("pipeline/output/output_20250519_142700_machine_learning_human/" + name + "/")

#valid_go_term_file = "pipeline/output/output_20250515_105213/Homo_sapiens/cellular_child_terms.txt"
valid_go_term_file = os.path.join(working_dir, "legend_content.txt")

# Parse the valid GO terms from the file
with open(valid_go_term_file, "r") as file:
    valid_go_terms = set(line.strip() for line in file)

# Read the feature matrix from the CSV file
feature_matrix_path = working_dir + "/feature_matrix_with_go_terms.csv"
df = pd.read_csv(feature_matrix_path, index_col=0)

# Drop the "Sequence" column if it exists
if "Sequence" in df.columns:
    df = df.drop(columns=["Sequence"])

# Create output directory for results
general_output_dir = os.path.join(working_dir, "go_term_rf_results10")
os.makedirs(general_output_dir, exist_ok=True)

go_ids = [
        #Go terms for: ER, Golgi, Ribosome, Mitochondria, Nucleus, Lysosome, Cell membrane, Cytoplasm
        "GO:0005783",
        "GO:0005794",
        "GO:0005840",
        "GO:0005739_no_cleavable_mts",
        "GO:0005739_cleavable_mts"
        "GO:0005634",
        "GO:0005764",
        "GO:0005886",
        "GO:0005737"
    ]
valid_go_terms = go_ids

for go_term in valid_go_terms:
    output_dir = os.path.join(general_output_dir, sanitize_filename(get_go_aspect(go_term, go_dag)))
    os.makedirs(output_dir, exist_ok=True)
    # Filter for current GO term (binary classification: this term vs. all others)
    df_filtered = df.copy()
    df_filtered["GO_Term_Binary"] = (df_filtered["GO_Term"] == go_term).astype(int)
    # Skip if not enough samples for both classes
    if df_filtered["GO_Term_Binary"].sum() < 5 or (df_filtered["GO_Term_Binary"] == 0).sum() < 5:
        continue

    X = df_filtered.drop(["GO_Term", "GO_Term_Binary"], axis=1)
    y = df_filtered["GO_Term_Binary"]

    # Count how many proteins have the current GO term
    proteins_with_go_term = df_filtered["GO_Term_Binary"].sum()

    # Add intercept for statsmodels
    X = drop_one_hot_collinear(X)
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
    rejected, pvals_corrected, _, _ = multipletests(feature_importance_df['Significance'], alpha=0.1, method='fdr_bh')
    feature_importance_df['FDR'] = pvals_corrected
    feature_importance_df['FDR_significant'] = rejected

    # Save features with their significance and FDR to file
    feature_importance_df.to_csv(os.path.join(output_dir, f"{go_term}_logreg_coefficients.csv"), index=False)

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

    plt.figure(figsize=(10, 6))
    sns.barplot(x='Coefficient', y='Feature', data=feature_importance_df.head(10))
    plt.title(f"Top 10 Logistic Regression Coefficients for GO Term: {go_term}")
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"{go_term}_logreg_coefficients_top10.png"))
    plt.close()

    # plot a volcano plot
    plt.figure(figsize=(10, 6))
    sns.scatterplot(x='Coefficient', y=-np.log10(feature_importance_df['Significance']), data=feature_importance_df, color='grey')
    plt.title(f"Volcano Plot for GO Term: {go_term}, {get_go_aspect(go_term, go_dag)}, n = {proteins_with_go_term}")
    plt.xlabel("Coefficient")
    plt.ylabel("-log10(Significance)")
    plt.axhline(y=-np.log10(0.05), color='r', linestyle='--')
    plt.axvline(x=0, color='g', linestyle='--')
    # Highlight significant features
    significant_features = feature_importance_df[feature_importance_df['FDR'] < 0.1]
    plt.scatter(significant_features['Coefficient'], -np.log10(significant_features['Significance']), color='red', label='Significant')

    # Use adjustText for non-overlapping labels
    texts = []
    for _, row in significant_features.iterrows():
        texts.append(
            plt.text(row['Coefficient'], -np.log10(row['Significance']), row['Feature'], fontsize=8, color='blue')
        )
    adjust_text(texts, arrowprops=dict(arrowstyle='->', color='black', lw=0.5))

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"{go_term}_volcano_plot.png"))
    plt.close()
    

if __name__ == "__main__":
# delete empty directories
    for root, dirs, files in os.walk(general_output_dir, topdown=False):
        for name in dirs:
            dir_path = os.path.join(root, name)
            if not os.listdir(dir_path):
                os.rmdir(dir_path)
                print(f"Deleted empty directory: {dir_path}")

