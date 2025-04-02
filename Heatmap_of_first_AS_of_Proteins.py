import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio import SeqIO
from scipy.stats import hypergeom

def fasta_to_dataframe(fasta_file):
    data = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        data.append([str(record.seq)])
    return pd.DataFrame(data, columns=['sequence'])


def count_instances_at_positions(array):
    rows, cols = array.shape
    print(rows, cols)
    dictlist = []
    for col in range(0,20): ######
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
amino_acid = np.array(["D", "E", "N", "Q", "Y", "H", "K", "R", "M", "L", "F", "I", "W", "S", "A", "T", "C", "P", "G", "V"])

input= ['output files/filtered_proteins_cleavable_mts.fasta', 
        'output files/filtered_proteins_no_cleavable_mts.fasta',
        'output files/filtered_proteins.fasta']

#wanted_result either "absolute" or "hgt"
wanted_result = "hgt"

all_arrays = [0, 0, 0]
all_counted_instances = [0, 0, 0]
amount_of_proteins = [0, 0, 0]
for file in range(len(input)):
    panda_df = fasta_to_dataframe(input[file])
    proteomes = list(panda_df['sequence'])
    #print(proteomes)

    proteome_array = np.zeros((len(proteomes), 20), dtype=object)
    for i, proteome in enumerate(proteomes):
        proteome = proteome[1:21] 
        for j in range(20):
            if j < len(proteome):
                proteome_array[i, j] = proteome[j]
            else:
                proteome_array[i, j] = 'X'
  #  print(proteome_array)

    counted_instances = count_instances_at_positions(proteome_array)
    rows, cols = proteome_array.shape
    amount_of_proteins[file] = rows
    all_counted_instances[file] = counted_instances
    #print(counted_instances)

if wanted_result == "absolute":
    for file in range(2):
        counted_instances = all_counted_instances[file]
        cols = len(position)
        rows = 19
        visual_array = np.zeros((cols, rows))
        for i in range(rows):
            for j in range(cols):
                if amino_acid[j] in counted_instances[i]:
                    visual_array[j, i] = counted_instances[i][amino_acid[j]]
        all_arrays[file]= visual_array
        #print(visual_array)

if wanted_result == "hgt":
    whole_set = all_counted_instances[2]
    for file in range(2):
        counted_instances = all_counted_instances[file]
   #     print(counted_instances)
        cols = len(position)
        rows = len(amino_acid)
        visual_array = np.zeros((cols, rows))
    #    print(visual_array)
        for i in range(rows):
            for j in range(cols):
                if amino_acid[j] in counted_instances[i]:
                    x = counted_instances[i][amino_acid[j]] #amount of amino acid in the subset
                    n = whole_set[i][amino_acid[j]] #amount of amino acid in the whole set
                    M = amount_of_proteins[2] #amount of all amino acids in the whole set
                    N = amount_of_proteins[file] #amount of all amino acids in the subset
                    p_value = hypergeom.pmf(x, M, n, N)
                    abs_log_val = abs(np.log10(p_value))
                    f_obs = x / N
                    f_exp = n / M
 #                   print(f_obs, f_exp)
                    result = 0
                    if f_obs < f_exp:
                        result = -1 * abs_log_val
                    else:
                        result = abs_log_val
                    visual_array[j, i] = result
 #       print(visual_array)
        all_arrays[file]= visual_array
        print(visual_array)
 #       print(visual_array.shape)


fig, ax = plt.subplots(ncols= 2)

ax[0].set_xticks(range(len(position)), labels=position, rotation=0, ha="right", rotation_mode="anchor")
ax[0].set_yticks(range(len(amino_acid)), labels=amino_acid, rotation=0, rotation_mode="anchor")
ax[0].set_yticklabels(amino_acid)
ax[0].set_xticklabels(position)

# Normalize the colormap across both graphs
vmin = min(np.min(all_arrays[0]), np.min(all_arrays[1]))
vmax = max(np.max(all_arrays[0]), np.max(all_arrays[1]))
cmap = "coolwarm"
pcm = ax[0].imshow(all_arrays[0], cmap=cmap, vmin=vmin, vmax=vmax)
pcm = ax[1].imshow(all_arrays[1], cmap=cmap, vmin=vmin, vmax=vmax)

ax[0].set_xlabel("Position")
ax[0].set_ylabel("Amino acids")
ax[0].set_title("Mitochondrial Proteins with MTS")

ax[1].set_xticks(range(len(position)), labels=position, rotation=0, ha="right", rotation_mode="anchor")
ax[1].set_yticks(range(len(amino_acid)), labels=amino_acid)
ax[1].set_xlabel("Position")
ax[1].set_ylabel("Amino acids")
ax[1].set_title("Mitochondrial Proteins without MTS")

'''
for i in range(len(position)):
    for j in range(len(amino_acid)):
        text = ax.text(j, i, visual_array[i, j], ha="center", va="center", color="black")'
'''
fig.tight_layout()
fig.colorbar(pcm, ax=ax, shrink=0.6)
plt.show()