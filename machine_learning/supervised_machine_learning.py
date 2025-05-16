import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelEncoder
import matplotlib.pyplot as plt
import seaborn as sns
import os
from goatools.obo_parser import GODag
from goatools.base import download_go_basic_obo



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
#working_dir = os.path.dirname("pipeline/output/output_20250514_134354/" + name + "/")
working_dir = os.path.dirname("pipeline/output/output_20250515_105213/" + name + "/")
valid_go_term_file = "pipeline/output/output_20250515_105213/Homo_sapiens/cellular_child_terms.txt"
# Parse the valid GO terms from the file
with open(valid_go_term_file, "r") as file:
    valid_go_terms = set(line.strip() for line in file)

feature_matrix_path = working_dir + "/feature_matrix_with_go_terms.csv"
# Read the feature matrix from the CSV file
df = pd.read_csv(feature_matrix_path, index_col=0)
# Drop the "Sequence" column if it exists
if "Sequence" in df.columns:
    df = df.drop(columns=["Sequence"])
print(df.head(10))

# Create output directory for results
general_output_dir = os.path.join(working_dir, "go_term_rf_results")
os.makedirs(general_output_dir, exist_ok=True)

for go_term in valid_go_terms:
    output_dir = os.path.join(general_output_dir, go_term)
    os. makedirs(output_dir, exist_ok=True)
    # Filter for current GO term (binary classification: this term vs. all others)
    df_filtered = df.copy()
    df_filtered["GO_Term_Binary"] = (df_filtered["GO_Term"] == go_term).astype(int)
    # Skip if not enough samples for both classes
    if df_filtered["GO_Term_Binary"].sum() < 5 or (df_filtered["GO_Term_Binary"] == 0).sum() < 5:
        continue

    X = df_filtered.drop(["GO_Term", "GO_Term_Binary"], axis=1)
    y = df_filtered["GO_Term_Binary"]

    # Encode labels if necessary
    if y.dtype == 'object':
        y = LabelEncoder().fit_transform(y)

    # Train/test split
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

    # Train Random Forest
    model = RandomForestClassifier(n_estimators=100, random_state=42)
    model.fit(X_train, y_train)

    # Feature importances
    importances = model.feature_importances_
    feature_importance_df = pd.DataFrame({
        'Feature': X.columns,
        'Importance': importances
    }).sort_values(by='Importance', ascending=False)

    # Save top feature to file
    top_feature = feature_importance_df.iloc[0]
    with open(os.path.join(output_dir, f"{get_go_aspect(go_term)}_top_feature.txt"), "w") as f:
        f.write(f"Top Feature: {top_feature['Feature']}\nImportance: {top_feature['Importance']}\n")

    # Plot and save feature importances
    plt.figure(figsize=(10, 6))
    sns.barplot(x='Importance', y='Feature', data=feature_importance_df)
    plt.title(f"Feature Importance for GO Term: {go_term}")
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"{go_term}_feature_importance.png"))
    plt.close()
