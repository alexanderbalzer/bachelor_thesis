import pandas as pd
from scipy.stats import ttest_ind
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
import os

def format_species_name(name: str) -> str:
    # Teile den Namen anhand des Unterstrichs
    parts = name.split("_")
    if name == "Homo_sapiens_isoforms":
        return "H. Sapiens with isoforms"
    if len(parts) != 2:
        raise ValueError("Name muss genau ein Unterstrich enthalten (Gattung_Art)")
    genus, species = parts
    # Kürze den Gattungsnamen auf den ersten Buchstaben + Punkt
    short_genus = genus[0] + "."
    
    # Setze alles in kursiv (z. B. für Markdown oder HTML)
    formatted = f"{short_genus} {species}"
    return formatted

eisenberg_hydrophobicity = {
    'A':  0.62,  # Alanine
    'R': -2.53,  # Arginine
    'N': -0.78,  # Asparagine
    'D': -0.90,  # Aspartic acid
    'C':  0.29,  # Cysteine
    'Q': -0.85,  # Glutamine
    'E': -0.74,  # Glutamic acid
    'G':  0.48,  # Glycine
    'H': -0.40,  # Histidine
    'I':  1.38,  # Isoleucine
    'L':  1.06,  # Leucine
    'K': -1.50,  # Lysine
    'M':  0.64,  # Methionine
    'F':  1.19,  # Phenylalanine
    'P':  0.12,  # Proline
    'S': -0.18,  # Serine
    'T': -0.05,  # Threonine
    'W':  0.81,  # Tryptophan
    'Y':  0.26,  # Tyrosine
    'V':  1.08   # Valine
}


fauchere_pliskat = {
    'A':  0.31, 'R': -1.01, 'N': -0.60,
    'D': -0.77, 'C':  1.54, 'Q': -0.22,
    'E': -0.64, 'G':  0.00, 'H':  0.13,
    'I':  1.80, 'L':  1.70, 'K': -0.99,
    'M':  1.23, 'F':  1.79, 'P':  0.72,
    'S': -0.04, 'T':  0.26, 'W':  2.25,
    'Y':  0.96, 'V':  1.22}
kyte_doolittle = {
    'A':  1.8,   # Alanine
    'R': -4.5,   # Arginine
    'N': -3.5,   # Asparagine
    'D': -3.5,   # Aspartate
    'C':  2.5,   # Cysteine
    'Q': -3.5,   # Glutamine
    'E': -3.5,   # Glutamate
    'G': -0.4,   # Glycine
    'H': -3.2,   # Histidine
    'I':  4.5,   # Isoleucine
    'L':  3.8,   # Leucine
    'K': -3.9,   # Lysine
    'M':  1.9,   # Methionine
    'F':  2.8,   # Phenylalanine
    'P': -1.6,   # Proline
    'S': -0.8,   # Serine
    'T': -0.7,   # Threonine
    'W': -0.9,   # Tryptophan
    'Y': -1.3,   # Tyrosine
    'V':  4.2    # Valine
}
hopp_woods = {
    'A': -0.5,   # Alanine
    'R':  3.0,   # Arginine
    'N':  0.2,   # Asparagine
    'D':  3.0,   # Aspartate
    'C': -0.3,   # Cysteine
    'Q':  0.2,   # Glutamine
    'E':  3.0,   # Glutamate
    'G':  0.0,   # Glycine
    'H': -0.5,   # Histidine
    'I': -1.8,   # Isoleucine
    'L': -1.8,   # Leucine
    'K':  3.0,   # Lysine
    'M': -1.3,   # Methionine
    'F': -2.5,   # Phenylalanine
    'P':  0.3,   # Proline
    'S':  0.3,   # Serine
    'T': -0.4,   # Threonine
    'W': -3.4,   # Tryptophan
    'Y': -2.3,   # Tyrosine
    'V': -1.5    # Valine
}

amino_acid_volumes = {
    'A':  88.6,   # Alanine
    'R': 173.4,   # Arginine
    'N': 114.1,   # Asparagine
    'D': 111.1,   # Aspartate
    'C': 108.5,   # Cysteine
    'Q': 143.8,   # Glutamine
    'E': 138.4,   # Glutamate
    'G':  60.1,   # Glycine
    'H': 153.2,   # Histidine
    'I': 166.7,   # Isoleucine
    'L': 166.7,   # Leucine
    'K': 168.6,   # Lysine
    'M': 162.9,   # Methionine
    'F': 189.9,   # Phenylalanine
    'P': 112.7,   # Proline
    'S':  89.0,   # Serine
    'T': 116.1,   # Threonine
    'W': 227.8,   # Tryptophan
    'Y': 193.6,   # Tyrosine
    'V': 140.0    # Valine
}
nat_specificity = {
    'A': 'NatA',   # Alanine
    'C': 'NatA',   # Cysteine
    'G': 'NatA',   # Glycine
    'S': 'NatA',   # Serine
    'T': 'NatA',   # Threonine
    'V': 'NatA',   # Valine
    'L': 'NatC',   # Leucine
    'I': 'NatC',   # Isoleucine
    'F': 'NatC',   # Phenylalanine
    'W': 'NatC',   # Tryptophan
    'Y': 'NatC',   # Tyrosine
    'M': 'Not acetylated',   # Methionine (retained, when P2 prevents iMet cleavage)
    'D': 'NatB',   # Aspartate
    'E': 'NatB',   # Glutamate
    'N': 'NatB',   # Asparagine
    'Q': 'NatB',   # Glutamine
    'K': 'Not acetylated',   # Lysine
    'R': 'Not acetylated',   # Arginine (borderline)
    'H': 'Not acetylated',   # Histidine
    'P': 'Not acetylated',  # Proline blocks NAT recognition (steric hindrance)
}
custom_values = {
    'D': -14.7637807665151,
    'E': -15.3294025001901,
    'N': -3.07926236466163,
    'Q': 1.03222407213378,
    'Y': -0.364775769848083,
    'H': -1.2449361016032,
    'K': -5.41036175294662,
    'R': -0.334479620065325,
    'M': -0.793970456423026,
    'L': 27.532794670514,
    'F': 3.45679037323112,
    'I': 1.43323617526077,
    'W': 8.84816338528988,
    'S': -5.19677910820452,
    'A': 25.1059608676786,
    'T': -1.38184389231379,
    'C': -1.42172035411789,
    'P': -2.45576810438212,
    'G': -8.22653087338568,
    'V': -2.12461123525388
}
helix_propensity_Chou_Fasman = {
    'A': 1.41,  # Alanine
    'R': 0.98,  # Arginine
    'N': 0.67,  # Asparagine
    'D': 0.90,  # Aspartate
    'C': 0.70,  # Cysteine
    'Q': 1.11,  # Glutamine
    'E': 1.51,  # Glutamate
    'G': 0.57,  # Glycine
    'H': 1.00,  # Histidine
    'I': 1.08,  # Isoleucine
    'L': 1.34,  # Leucine
    'K': 1.07,  # Lysine
    'M': 1.45,  # Methionine
    'F': 1.12,  # Phenylalanine
    'P': 0.59,  # Proline (helix breaker)
    'S': 0.77,  # Serine
    'T': 0.83,  # Threonine
    'W': 1.14,  # Tryptophan
    'Y': 0.61,  # Tyrosine
    'V': 0.91   # Valine
}

amino_acids = list(custom_values.keys())

organism = 'Saccharomyces_cerevisiae'  # Change this to the desired organism
#organism = 'Homo_sapiens'  # Change this to the desired organism

# Load the feature matrix (replace 'feature_matrix.csv' with your file path)
file_path = f'pipeline/output/output_20250617_183139_latest_ML/{organism}/feature_matrix_with_go_terms.csv'
#file_path = 'pipeline/output/output_20250617_183139_latest_ML/Saccharomyces_cerevisiae/feature_matrix_with_go_terms.csv'
data = pd.read_csv(file_path)
cyto_nuclear_df = data.copy()

# Filter for mitochondrial proteins, take their sequence and map if they fit to a specific mask

mito_df = data[(data['GO_Term'] == 'GO:0005739') & (data['Hydrophobic Moment'] >= 0.5)]
cyto_nuclear_df = cyto_nuclear_df[(cyto_nuclear_df['GO_Term'] == 'cyto_nuclear')]
# Categorize by the second amino acid
categories1 = pd.DataFrame({
    'aas': cyto_nuclear_df['Sequence'].str[1:3].fillna('Unknown')
})
categories2 = pd.DataFrame({
    'aas': mito_df['Sequence'].str[1:3].fillna('Unknown')
})
nat = { 
    "MA": 'NatA',
    "MS": 'NatA',
    "MT": 'NatA',
    "MG": 'NatA',
    "MV": 'NatA',
    "MC": 'NatA',
    "ML": 'NatC',
    "MF": 'NatC',
    "MI": 'NatC',
    "MW": 'NatC',
    "MY": 'NatC',
    "MK": 'NatC',
    "MD": 'NatB',
    "ME": 'NatB',
    "MN": 'NatB',
    "MQ": 'NatB'
}


hydrophobic = set('MALFIWVC')
neutral = set('STGPYH')
hydrophilic = set('NQEDRK')
def get_category(seq):
    if len(seq) < 2:
        return 'Unknown'
    masked_seq = ''
    for aa in seq:
        if aa in hydrophobic:
            masked_seq += 'H'
        elif aa in neutral:
            masked_seq += 'N'
        elif aa in hydrophilic:
            masked_seq += 'P'
        else:
            masked_seq += 'X'
    return masked_seq

categories1['Category'] = categories1['aas'].apply(get_category)
categories2['Category'] = categories2['aas'].apply(get_category)

categories1['Categories_grouped'] = categories1['Category']
categories2['Categories_grouped'] = categories2['Category']

# Count occurrences of each category
value_counts = categories1['Categories_grouped'].value_counts()
value_counts2 = categories2['Categories_grouped'].value_counts()
# Create subplots for the pie charts
fig, axes = plt.subplots(1, 2, figsize=(16, 8))

# Pie chart for the first value counts
axes[0].pie(value_counts, labels=value_counts.index, autopct='%1.1f%%', startangle=90)
axes[0].set_title(f'Distribution of Second aa in cyto_nuclear (n = {len(cyto_nuclear_df)})')

# Pie chart for the second value counts
axes[1].pie(value_counts2, labels=value_counts2.index, autopct='%1.1f%%', startangle=90)
axes[1].set_title(f'Distribution of Second aa in mito (n = {len(mito_df)})')

# Adjust layout and show the figure
plt.tight_layout()

working_dir = 'pipeline/output/output_20250617_183904'
output_dir = os.path.join(working_dir, 'plots')
plt.savefig(os.path.join(output_dir, f'distribution_of_second_aa_in_cyto_nuclear_and_mito{organism}.pdf'), dpi=300)
plt.show()
exit()

# Map the hydrophobicity values to the amino acids in value_counts2
# Prepare data for 3D scatter plot
organism_names = [
    "Homo_sapiens","Mus_musculus", "Rattus_norvegicus", "Danio_rerio",
    "Caenorhabditis_elegans", "Drosophila_melanogaster", "Arabidopsis_thaliana", 
    "Saccharomyces_cerevisiae"]
for organism in organism_names:
    working_dir = 'pipeline/output/output_20250617_183904'
    output_dir = os.path.join(working_dir, 'plots')
    os.makedirs(output_dir, exist_ok=True)
    visual_array_path = os.path.join(working_dir, 'visual_array_organisms.csv')
    hgt_array_df = pd.read_csv(visual_array_path, index_col=0)
    hgt_values = hgt_array_df.loc[organism]
    hydrophobicity_values = [eisenberg_hydrophobicity.get(aa, np.nan) for aa in amino_acids]
    helix_propensity_values = [helix_propensity_Chou_Fasman.get(aa, np.nan) for aa in amino_acids]
    nat_substrate_types = [nat_specificity.get(aa, 'Unknown') for aa in amino_acids]
    custom_values_list = [hgt_values[aa] if aa in hgt_values else np.nan for aa in amino_acids]


    # Create a DataFrame for the 3D scatter plot
    scatter_data = pd.DataFrame({
        'Amino Acid': amino_acids,
        'Hydrophobicity (Eisenberg)': hydrophobicity_values,
        'Helix Propensity (Chou-Fasman)': helix_propensity_values,
        'HGT score': custom_values_list,
        'NAT Substrate Type': nat_substrate_types
    })
    # 


    fig2 = px.scatter(
        scatter_data,
        x='Hydrophobicity (Eisenberg)',
        y='Helix Propensity (Chou-Fasman)',
        color='HGT score',
        color_continuous_scale='RdBu_r',
        range_color=[-20, 20],
        symbol='NAT Substrate Type',
        text='Amino Acid',
        title=f'Amino Acid Properties, HGT Scores in proteins with MTS in {format_species_name(organism)}'
    )
    fig2.update_traces(
        marker=dict(
            size=16,
            line=dict(color='black', width=1)
        ),
        textposition='middle left'
    )
    fig2.update_layout(
        legend=dict(
            x=0,
            y=1,
            xanchor='left',
            yanchor='top'
        ),
        coloraxis_colorbar=dict(
            tickvals=[-20, 0, 20],
            ticktext=['-20', '0', '20'],
            title='HGT',
            len=0.5
        )
    )

    # Remove the original legend for NAT substrate type
    #fig2.for_each_trace(lambda t: t.update(showlegend=False))

    # Add a dummy trace for each NAT substrate type with white marker and black outline
    '''symbols = {
        'NatA': 'circle',
        'NatB': 'square',
        'NatC': 'diamond',
        'Not acetylated': 'x',
        'Unknown': 'triangle-up'
    }
    for nat_type in scatter_data['NAT Substrate Type'].unique():
        fig2.add_trace(go.Scatter(
            x=[None], y=[None],
            mode='markers',
            marker=dict(
                symbol=symbols.get(nat_type, 'circle'),
                color='white',
                size=16,
                line=dict(color='black', width=1)
            ),
            name=nat_type,
            showlegend=True
        ))'''
    output_file_html = os.path.join(output_dir, f"hydrophobicity_vs_helix_propensity.html")
    output_file_png = os.path.join(output_dir, f"hydrophobicity_vs_helix_propensity{organism}.pdf")
    
    # Save the figure as a PDF
    fig2.write_image(output_file_png, scale=5)

'''
# Normalize hydrophobicity and helix propensity
scatter_data['Normalized Hydrophobicity'] = (scatter_data['Hydrophobicity'] - scatter_data['Hydrophobicity'].min()) / (scatter_data['Hydrophobicity'].max() - scatter_data['Hydrophobicity'].min())
scatter_data['Normalized Helix Propensity'] = (scatter_data['Helix Propensity'] - scatter_data['Helix Propensity'].min()) / (scatter_data['Helix Propensity'].max() - scatter_data['Helix Propensity'].min())

# Create a new value by adding normalized hydrophobicity and helix propensity
scatter_data['Combined Score'] = scatter_data['Normalized Hydrophobicity'] + scatter_data['Normalized Helix Propensity']
# Plot HGT score against the new combined score
'''


# Create subplots for the pie charts
fig, axes = plt.subplots(1, 2, figsize=(16, 8))

# Pie chart for the first value counts
axes[0].pie(value_counts, labels=value_counts.index, autopct='%1.1f%%', startangle=90)
axes[0].set_title('Distribution of Second aa in cyto_nuclear')

# Pie chart for the second value counts
axes[1].pie(value_counts2, labels=value_counts2.index, autopct='%1.1f%%', startangle=90)
axes[1].set_title('Distribution of Second aa in mito')

# Adjust layout and show the figure
plt.tight_layout()
# Perform a t-test between the two distributions
# Ensure both distributions have the same index to align them
common_categories = value_counts.index.intersection(value_counts2.index)
dist1 = value_counts[common_categories]
dist2 = value_counts2[common_categories]

# Perform the t-test
t_stat, p_value = ttest_ind(dist1, dist2, equal_var=False)

# Print the results
print(f"T-statistic: {t_stat}")
print(f"P-value: {p_value}")

# Count occurrences of each unique value in the 'H_first_two_aa_neo_N' column
H_mito = mito_df['H_first_two_aa_neo_N']
H_prot = cyto_nuclear_df['H_first_two_aa_neo_N']

# Drop NaN values for plotting
H_mito_clean = H_mito.dropna()
H_prot_clean = H_prot.dropna()

# Create a DataFrame for violin plot
violin_data = pd.DataFrame({
    'Category': ['H_mito'] * len(H_mito_clean) + ['H_prot'] * len(H_prot_clean),
    'Values': pd.concat([H_mito_clean, H_prot_clean])
})

# Plot violin plot
plt.figure(figsize=(10, 6))
sns.violinplot(x='Category', y='Values', data=violin_data)
plt.title('Hydrophobicity of First Two Amino Acids of N-terminus')
plt.ylabel('Hydrophobicity')
plt.tight_layout()
# Perform a t-test for H_mito and H_prot
# Drop NaN values to ensure valid comparison

# Perform the t-test
t_stat_H, p_value_H = ttest_ind(H_mito_clean, H_prot_clean, equal_var=False)

# Print the results
print(f"T-statistic for H_mito vs H_prot: {t_stat_H}")
print(f"P-value for H_mito vs H_prot: {p_value_H}")
