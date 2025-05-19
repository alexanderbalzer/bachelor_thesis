import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
from goatools.obo_parser import GODag
from goatools.base import download_go_basic_obo
import statsmodels.api as sm
from sklearn.preprocessing import StandardScaler, LabelEncoder

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

# Read the feature matrix from the CSV file
feature_matrix_path = working_dir + "/feature_matrix_with_go_terms.csv"
df = pd.read_csv(feature_matrix_path, index_col=0)

# Drop the "Sequence" column if it exists
if "Sequence" in df.columns:
    df = df.drop(columns=["Sequence"])

# Encode GO_Term as integer labels
le = LabelEncoder()
df["GO_Term_Label"] = le.fit_transform(df["GO_Term"])

X = df.drop(columns=["GO_Term", "GO_Term_Label"])
y = df["GO_Term_Label"]

# Z-transformation (standardization)
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)
X_scaled = pd.DataFrame(X_scaled, columns=X.columns, index=X.index)
X_scaled = sm.add_constant(X_scaled)

# Fit multinomial logistic regression
model = sm.MNLogit(y, X_scaled)
result = model.fit(method="newton", maxiter=100, disp=0)
# Create output directory for results
output_dir = os.path.join(working_dir, "go_term_multiclass_logreg_results")
os.makedirs(output_dir, exist_ok=True)

# Get coefficients (result.params is shape [n_classes-1, n_features+1])
for class_idx, go_term_label in enumerate(le.classes_):
    if class_idx == 0:
        continue  # Reference class, no coefficients
    coefs = result.params.iloc[class_idx - 1].drop("const", errors="ignore")
    print(f"GO term: {go_term_label}, class_idx: {class_idx}, coefs: {coefs.head()}")  # Debug print
    abs_coefs = np.abs(coefs)
    top10_idx = abs_coefs.sort_values(ascending=False).index[:10]
    top10_features = coefs[top10_idx]
    print(f"Top 10 features for {go_term_label}: {top10_features}")  # Debug print

    # Save top 10 features to file
    with open(os.path.join(output_dir, f"{go_term_label}_top_10_features.txt"), "w") as f:
        for feat, coef in top10_features.items():
            f.write(f"{feat}\tCoefficient: {coef}\n")

    # Plot and save
    plt.figure(figsize=(10, 6))
    sns.barplot(x=top10_features.values, y=top10_features.index)
    plt.title(f"Top 10 Logistic Regression Coefficients for GO Term: {go_term_label}")
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"{go_term_label}_logreg_coefficients_top10.png"))
    plt.close()

