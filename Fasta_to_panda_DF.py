import fastaparser as fp
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import pandas as pd


def fasta_to_dataframe(fasta_file):
    data = []
    with open(fasta_file, 'r') as file:
        for sequence in fp.Reader(file):
            data.append([sequence.id, sequence.sequence])
        return pd.DataFrame(data, columns=['id', 'sequence'])

def count_instances_at_positions(array):
    print(array.shape)
    for col in range(cols) < 20:
        count_dict = {}
        dictlist = []
        for row in range(rows):
            value = array[row, col]
            if (col, value) in count_dict:
                count_dict[(col, value)] += 1
            else:
                count_dict[(col, value)] = 1
        dictlist.append[count_dict]
    return dictlist


position = np.array(["2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20"])
amino_acid = np.array(["D", "E", "N", "Q", "Y", "H", "K", "R", "M", "L", "F", "I", "W", "S", "T", "C", "P", "G", "V"])

panda_df = fasta_to_dataframe('P53_HUMAN.fasta')
for i in range(len(panda_df)):
    proteomes = ['']*len(panda_df)
    proteomes[i] = [panda_df.iat[i,1]]
#print(proteomes)
counted_instances = count_instances_at_positions(proteomes)
proteome_array = np.array(counted_instances)

print(proteome_array)
fig, ax = plt.subplots()
im = ax.imshow(proteome_array)
ax.set_xticks(range(len(position)), labels=position, rotation=45, ha="right", rotation_mode="anchor")
ax.set_yticks(range(len(amino_acid)), labels=amino_acid)

for i in range(len(position)):
    for j in range(len(amino_acid)):
        text = ax.text(j, i, proteome[i, j], ha="center", va="center", color="w")
ax.set_title("Mitochondrial Proteins with MTS")
fig.tight_layout()
plt.show()