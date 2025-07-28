import os
import pandas as pd

name = "Homo_sapiens"  # Example organism name

# Set the working directory
working_dir = os.path.dirname("pipeline/output/output_20250519_142700_machine_learning_human/" + name + "/")

# Read the feature matrix from the CSV file
feature_matrix_path = working_dir + "/feature_matrix_with_go_terms.csv"
df = pd.read_csv(feature_matrix_path)
uniprot_ids = ["TOM5", "TOM6", "TOM7",
"TOM20", "TOM22", "TOM40",
"NACA", "BTF3", "PINK1", "VDAC2"]
uniprot_ids = [prot_id + "_HUMAN" for prot_id in uniprot_ids]
print(uniprot_ids)
# Filter the dataframe for the specified proteins
# Ensure the column name matches the actual column name in the CSV file
filtered_df = df[df['protein_id_human'].isin(uniprot_ids)]

# Create a new dataframe with the sequences split at the mpp_cleavage_pos
split_sequences = []
for _, row in filtered_df.iterrows():
    cleavage_pos = int(row['MPP_cleavage_position'])
    sequence = row['Sequence']
    split_sequences.append({
        'protein_id': row['protein_id_human'],
        'presequence': sequence[:cleavage_pos],
        'mature_sequence': sequence[cleavage_pos:]
    })

# Convert the split sequences into a new dataframe
split_sequences_df = pd.DataFrame(split_sequences)

# Display or save the resulting dataframe
print(split_sequences_df)
