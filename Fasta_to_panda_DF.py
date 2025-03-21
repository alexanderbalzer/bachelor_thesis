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
    print(dictlist)
    return dictlist


position = np.array(["2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20"])
amino_acid = np.array(["D", "E", "N", "Q", "Y", "H", "K", "R", "M", "L", "F", "I", "W", "S", "T", "C", "P", "G", "V"])

panda_df = fasta_to_dataframe('P53_HUMAN.fasta')
proteomes = list(panda_df['sequence'])
#print(proteomes)

proteome_array = np.array([list(seq[:20]) for seq in proteomes])
print(proteome_array)

counted_instances = count_instances_at_positions(proteome_array)
print(counted_instances)


#fig, ax = plt.subplots()
#im = ax.imshow(counted_instances, cmap="viridis")
#ax.set_xticks(range(len(position)), labels=position, rotation=45, ha="right", rotation_mode="anchor")
#ax.set_yticks(range(len(amino_acid)), labels=amino_acid)

#for i in range(len(position)):
#    for j in range(len(amino_acid)):
#        text = ax.text(j, i, proteome[i, j], ha="center", va="center", color="w")
#ax.set_title("Mitochondrial Proteins with MTS")
#fig.tight_layout()
#plt.show()