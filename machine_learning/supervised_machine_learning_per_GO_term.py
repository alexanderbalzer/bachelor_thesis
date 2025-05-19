import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
from goatools.obo_parser import GODag
from goatools.base import download_go_basic_obo
import statsmodels.api as sm
from sklearn.preprocessing import StandardScaler

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
        return "Unknown"
    

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
print(df.head(10))

# Create output directory for results
general_output_dir = os.path.join(working_dir, "go_term_rf_results2")
os.makedirs(general_output_dir, exist_ok=True)

for go_term in valid_go_terms:
    output_dir = os.path.join(general_output_dir, get_go_aspect(go_term, go_dag))
    os.makedirs(output_dir, exist_ok=True)
    # Filter for current GO term (binary classification: this term vs. all others)
    df_filtered = df.copy()
    df_filtered["GO_Term_Binary"] = (df_filtered["GO_Term"] == go_term).astype(int)
    # Skip if not enough samples for both classes
    if df_filtered["GO_Term_Binary"].sum() < 5 or (df_filtered["GO_Term_Binary"] == 0).sum() < 5:
        continue

    X = df_filtered.drop(["GO_Term", "GO_Term_Binary"], axis=1)
    y = df_filtered["GO_Term_Binary"]

    # Add intercept for statsmodels
    X = sm.add_constant(X)

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
    go_term = get_go_aspect(go_term, go_dag)
    # Get coefficients (excluding intercept)
    coefs = result.params.drop("const")
    feature_importance_df = pd.DataFrame({
        'Feature': coefs.index,
        'Coefficient': coefs.values,
        'AbsCoefficient': np.abs(coefs.values)
    }).sort_values(by='AbsCoefficient', ascending=False)

    # Save top feature to file
    top_feature = feature_importance_df.iloc[0]
    with open(os.path.join(output_dir, f"{go_term}_top_feature.txt"), "w") as f:
        f.write(f"Top Feature: {top_feature['Feature']}\nCoefficient: {top_feature['Coefficient']}\n")

    # Get top 10 features by absolute coefficient
    top_10_features = feature_importance_df.head(10)

    # Save top 10 features to file
    with open(os.path.join(output_dir, f"{go_term}_top_10_features.txt"), "w") as f:
        for idx, row in top_10_features.iterrows():
            f.write(f"{row['Feature']}\tCoefficient: {row['Coefficient']}\n")

    # Plot and save feature importances
    '''plt.figure(figsize=(10, 6))
    sns.barplot(x='Coefficient', y='Feature', data=feature_importance_df)
    plt.title(f"Logistic Regression Coefficients for GO Term: {go_term}")
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"{go_term}_logreg_coefficients.png"))
    plt.close()'''

    plt.figure(figsize=(10, 6))
    sns.barplot(x='Coefficient', y='Feature', data=top_10_features)
    plt.title(f"Top 10 Logistic Regression Coefficients for GO Term: {go_term}")
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"{go_term}_logreg_coefficients_top10.png"))
    plt.close()
    
