import pandas as pd
import matplotlib.pyplot as plt
from mygene import MyGeneInfo
from gprofiler import GProfiler
from tqdm import tqdm
from Bio import Entrez
import requests
import mygene

def get_uniprot_disease_annotation(uniprot_id):
    """
    Get disease-related information from UniProt for a given UniProt ID.
    """
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"
    try:
        response = requests.get(url)
        if response.status_code != 200:
            return f"Error fetching UniProt data: {response.status_code}"

        data = response.json()

        diseases = []
        for comment in data.get("comments", []):
            if comment["commentType"] == "DISEASE":
                disease_info = {
                    "disease": comment.get("disease", {}).get("diseaseId"),
                    "name": comment.get("disease", {}).get("diseaseAccession", ""),
                    "description": comment.get("disease", {}).get("description"),
                    "evidence": [ev["evidenceCode"] for ev in comment.get("evidences", [])]
                }
                diseases.append(disease_info)

        if not diseases:
            return "No disease annotation found."

        return diseases

    except Exception as e:
        return f"Error: {str(e)}"

def get_uniprot_function(uniprot_id):
    """
    Retrieve the functional description of a protein from UniProt using its UniProt ID.
    """
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"
    try:
        response = requests.get(url)
        if response.status_code != 200:
            return f"Error: Failed to fetch data for {uniprot_id} (HTTP {response.status_code})"

        data = response.json()
        function_comments = [
            comment for comment in data.get("comments", [])
            if comment.get("commentType") == "FUNCTION"
        ]

        if not function_comments:
            return "No functional annotation found."

        # Some entries may have multiple texts
        functions = []
        for entry in function_comments:
            texts = entry.get("texts", [])
            for text in texts:
                functions.append(text.get("value"))

        return functions
    except Exception as e:
        return f"Error retrieving function: {str(e)}"

# Specify the path to your CSV file
csv_file_path = 'pipeline/output/output_20250519_142700_machine_learning_human/Homo_sapiens/feature_matrix_with_go_terms.csv'
output_file = 'pipeline/output/output_20250519_142700_machine_learning_human/Homo_sapiens/'

# Read the CSV file into a pandas DataFrame
df = pd.read_csv(csv_file_path)

# Specify the column you want to process
diff_if_not_that_nat_substrate = df["diff_if_not_that_nat_substrate"]
diff_if_huntington = pd.Series([value if go_term != "cyto_nuclear" else None for value, go_term in zip(df["diff_if_huntington"], df['GO_Term'])])
sequence = df['Sequence']

# Extract the second amino acid from each sequence
second_aa = sequence.str[1]
# Filter proteins with positive diff_if_huntington
positive_diff_proteins = df[(df['GO_Term'] == "GO:0005739") & (df['GO_Term'] != 'Multiple') & (df['NAT_NatA/D'] == 1) & (df['Hydrophobic Moment'] > 0.3) & (df['Hydrophobic Moment'] -  df['diff_if_huntington'] < 0.3)]

for i, row in tqdm(positive_diff_proteins.iterrows(), total=len(positive_diff_proteins)):
    uniprot_id = row['protein_id']
    function = get_uniprot_function(uniprot_id)
    positive_diff_proteins.loc[i, 'Function'] = ",".join(function) if isinstance(function, list) else function
for i, row in tqdm(positive_diff_proteins.iterrows(), total=len(positive_diff_proteins)):
    uniprot_id = row['protein_id']  
    disease_annotations = get_uniprot_disease_annotation(uniprot_id)
    positive_diff_proteins.loc[i, 'Disease_Annotations'] = ",".join(
        [annotation['description'] or "" for annotation in disease_annotations
         if 'description' in annotation]
    )

# Save the filtered proteins to a new CSV file
positive_diff_proteins.to_csv(output_file + 'proteins_with_broken_mts_when_huntington.csv', index=False)

# For diff_if_not_that_nat_substrate
plt.figure(figsize=(10, 6))
plt.violinplot(
    [diff_if_not_that_nat_substrate[second_aa == aa].dropna().values for aa in second_aa.unique()],
    showmeans=True, showextrema=True
)
plt.xticks(ticks=range(1, len(second_aa.unique()) + 1), labels=second_aa.unique(), rotation=90)
plt.xlabel('Second Amino Acid')
plt.ylabel('Diff If Not That Nat Substrate')
plt.title('Diff If Not That Nat Substrate vs Second Amino Acid')
plt.tight_layout()

# For diff_if_huntington
plt.figure(figsize=(10, 6))
plt.violinplot(
    [diff_if_huntington[second_aa == aa].dropna().values for aa in second_aa.unique()],
    showmeans=True, showextrema=True
)
plt.xticks(ticks=range(1, len(second_aa.unique()) + 1), labels=second_aa.unique(), rotation=90)
plt.xlabel('Second Amino Acid')
plt.ylabel('diff_if_huntington')
plt.title('diff_if_huntington vs Second Amino Acid')
plt.tight_layout()
plt.show()