import pandas as pd
import os
# Load the CSV file
file_path = 'pipeline/output/output_20250519_142700_machine_learning_human/Homo_sapiens/feature_matrix_with_go_terms.csv'  # Replace with your CSV file path
data = pd.read_csv(file_path, index_col=0)

# Define the conditions for splitting
mito_protein = []
er_protein = []
cytosol_protein = []
somethings_wrong = []

for i, row in data.iterrows():
    if row['GO_Term'] == 'GO:0005739':
        mito_protein.append(row)
    elif row['GO_Term'] == 'GO:0005783':
        er_protein.append(row)
    elif row['GO_Term'] == 'cyto_nuclear':
        cytosol_protein.append(row)
    else:
        somethings_wrong.append(row)


# Convert the lists to DataFrames
df1 = pd.DataFrame(mito_protein)
df2 = pd.DataFrame(er_protein)
df3 = pd.DataFrame(cytosol_protein)
df4 = pd.DataFrame(somethings_wrong)
print(f"mito_proteins: {len(df1)}")
print(f"ER proteins: {len(df2)}")
print(f"cytosolic proteins: {len(df3)}")
print(f"proteins to er and mito: {len(df4)}")
print(f"sum: {len(df1)+ len(df2) + len(df3) + len(df4)}")
print(f"total: {len(data)}")

output_dir = 'pipeline/output/output_20250519_142700_machine_learning_human/Homo_sapiens/try'
os.makedirs(output_dir, exist_ok=True)
# Optionally, save the DataFrames to new CSV files
df1.to_csv(os.path.join(output_dir, "mito_proteins.csv"), index=False)
df2.to_csv(os.path.join(output_dir, "er_proteins.csv"), index=False)
df3.to_csv(os.path.join(output_dir, "cytosolic_proteins.csv"), index=False)