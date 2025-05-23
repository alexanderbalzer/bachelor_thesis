import pandas as pd
import matplotlib.pyplot as plt

# Specify the path to your CSV file
csv_file_path = 'pipeline/output/output_20250519_142700_machine_learning_human/Homo_sapiens/feature_matrix.csv'
output_file = 'pipeline/output/output_20250519_142700_machine_learning_human/Homo_sapiens/'

# Read the CSV file into a pandas DataFrame
df = pd.read_csv(csv_file_path)

# Specify the column you want to process
cleavage_pos = 'predicted_cleavage_position'
signalP = 'signalP_cleavage_probability'
MitoFates = 'cleavable_mts'

new_df = pd.DataFrame(columns=[cleavage_pos, signalP, MitoFates])
new_df = df[[cleavage_pos, signalP, MitoFates]].dropna()
print(new_df.head(10))
# Filter rows based on conditions and populate new_df
rows = []
for i, row in df.iterrows():
    if not pd.isna(row[cleavage_pos]) and not pd.isna(row[signalP]) and not pd.isna(row[MitoFates]):
        rows.append({
            cleavage_pos: row[cleavage_pos],
            signalP: row[signalP],
            MitoFates: row[MitoFates]
        })
new_df = pd.DataFrame(rows)

df_signalP = []
df_MitoFates = []
for i, row in new_df.iterrows():
    if row[cleavage_pos] == 24.016427:
        continue
    else:
        if row[signalP] > row[MitoFates]:
            if row[signalP] < 0.5:
                continue
            else:
                df_signalP.append(row[cleavage_pos])
        else:
            if row[MitoFates] < 0.5:
                continue
            else:
                df_MitoFates.append(row[cleavage_pos])

# Convert lists to DataFrames
df_signalP = pd.DataFrame(df_signalP, columns=[cleavage_pos])
df_MitoFates = pd.DataFrame(df_MitoFates, columns=[cleavage_pos])


df_signalP[cleavage_pos] = df_signalP[cleavage_pos].astype(float)
df_MitoFates[cleavage_pos] = df_MitoFates[cleavage_pos].astype(float)

# Plot histograms for the two DataFrames
plt.figure(figsize=(12, 6))

# Histogram for df_signalP
plt.subplot(1, 2, 1)
plt.hist(df_signalP[cleavage_pos], bins=30, color='blue', alpha=0.7)
plt.title('Histogram of SignalP Cleavage Positions')
plt.xlabel('Cleavage Position')
plt.ylabel('Frequency')

# Mark the mean for df_signalP
mean_signalP = df_signalP[cleavage_pos].mean()
plt.axvline(mean_signalP, color='red', linestyle='dashed', linewidth=1, label=f'Mean: {mean_signalP:.2f}')
plt.legend()

# Histogram for df_MitoFates
plt.subplot(1, 2, 2)
plt.hist(df_MitoFates[cleavage_pos], bins=30, color='green', alpha=0.7)
plt.title('Histogram of MitoFates Cleavage Positions')
plt.xlabel('Cleavage Position')
plt.ylabel('Frequency')

# Mark the mean for df_MitoFates
mean_MitoFates = df_MitoFates[cleavage_pos].mean()
plt.axvline(mean_MitoFates, color='red', linestyle='dashed', linewidth=1, label=f'Mean: {mean_MitoFates:.2f}')
plt.legend()

# Show the plots
plt.tight_layout()

# Save the histograms to files
plt.savefig(output_file + 'histograms_cleavage_pos_mts_ER.pdf', dpi=300)
plt.show()


