import pandas as pd
import matplotlib.pyplot as plt
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
        return go_id


# Load the GO DAG
go_dag = load_obo()

# Specify the path to your CSV file
csv_file_path = 'pipeline/output/output_20250519_142700_machine_learning_human/Homo_sapiens/feature_matrix_with_go_terms.csv'
output_file = 'pipeline/output/output_20250519_142700_machine_learning_human/Homo_sapiens/boxplots/'
os.makedirs(output_file, exist_ok=True)

# Read the CSV file into a pandas DataFrame
df = pd.read_csv(csv_file_path)

# Retrieve unique GO IDs from the DataFrame
go_ids = set()
for go_terms in df['GO_Term'].dropna():
    go_ids.update(go_terms.split(';'))

# Get all numeric columns (excluding index and GO_Term)
numeric_cols = df.drop(columns=["Sequence", "protein_id", "GO_Term"])
amount_of_features = numeric_cols.shape[1]
for i, feature in enumerate(numeric_cols):
    print(f"running feature {feature} ({i+1}/{amount_of_features})")
    data = []
    labels = []
    for go_id in go_ids:
        values = df[df['GO_Term'].str.contains(go_id, na=False)][feature].dropna()
        if len(values) > 0:
            data.append(values)
            labels.append(get_go_aspect(go_id, go_dag))
    if data:
        plt.figure(figsize=(max(8, len(labels)*1.2), 6))
        plt.boxplot(data, tick_labels=labels, patch_artist=True,
                    boxprops=dict(facecolor='lightblue'), medianprops=dict(color='red'))
        # Plot dots for each value (jittered)
        '''for i, values in enumerate(data):
            plt.scatter([i+1]*len(values), values, color='black', alpha=0.5, s=10, zorder=2)'''
        plt.title(f'{feature} Distribution by GO Term')
        plt.ylabel(feature)
        plt.xticks(rotation=45, ha='right')
        plt.grid(axis='y', linestyle='--', alpha=0.7)
        plt.tight_layout()
        plt.savefig(os.path.join(output_file, f'{feature}_boxplot_by_go_term.png'))
        plt.close()

