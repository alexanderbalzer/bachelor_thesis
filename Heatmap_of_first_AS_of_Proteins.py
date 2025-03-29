import fastaparser as fp
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import pandas as pd
from Bio import SeqIO

def fasta_to_dataframe(fasta_file):
    data = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        data.append([str(record.seq)])
    return pd.DataFrame(data, columns=['sequence'])


def count_instances_at_positions(array):
    rows, cols = array.shape
    print(rows, cols)
    dictlist = []
    for col in range(0,19):
        count_dict = {}
        for row in range(rows):
            value = array[row, col]
            if (value) in count_dict:
                count_dict[value] += 1
            else:
                count_dict[value] = 1
        dictlist.append(count_dict)
    #print(dictlist)
    return dictlist


position = np.array(["2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20"])
amino_acid = np.array(["D", "E", "N", "Q", "Y", "H", "K", "R", "M", "L", "F", "I", "W", "S", "T", "C", "P", "G", "V"])

panda_df = fasta_to_dataframe('output files/filtered_proteins_no_cleavable_mts.fasta')
proteomes = list(panda_df['sequence'])
#print(proteomes)

proteome_array = np.zeros((len(proteomes), 19), dtype=object)
for i, proteome in enumerate(proteomes):
    proteome = proteome[1:20]  # Extract amino acids from 2 to 20
    for j in range(19):
        if j < len(proteome):
            proteome_array[i, j] = proteome[j]
        else:
            proteome_array[i, j] = 'X'
#print(proteome_array)

counted_instances = count_instances_at_positions(proteome_array)
#print(counted_instances)

cols = len(amino_acid)
rows = 19
visual_array = np.zeros((rows, cols))
for i in range(rows):
    for j in range(cols):
        if amino_acid[j] in counted_instances[i]:
            visual_array[j, i] = counted_instances[i][amino_acid[j]]
#print(visual_array)

fig, ax = plt.subplots()
im = ax.imshow(visual_array, cmap="viridis")
ax.set_xticks(range(len(position)), labels=position, rotation=45, ha="right", rotation_mode="anchor")
ax.set_yticks(range(len(amino_acid)), labels=amino_acid)

for i in range(len(position)):
    for j in range(len(amino_acid)):
        text = ax.text(j, i, visual_array[i, j], ha="center", va="center", color="black")
ax.set_title("Mitochondrial Proteins without MTS")
fig.tight_layout()
plt.show()