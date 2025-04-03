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
    print("amount of proteins: ", rows)
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

Name = ["human"]
name = Name[0]

position = np.array(["2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20"])
amino_acid = np.array(["D", "E", "N", "Q", "Y", "H", "K", "R", "M", "L", "F", "I", "W", "S", "A", "T", "C", "P", "G", "V"])

#input= ["output files/filtered_proteins_cleavable_mts_for_" + str(name) + "_2.fasta", 
 #       "output files/filtered_proteins_no_cleavable_mts_for_" + str(name) + "_2.fasta",
  #      "output files/filtered_proteins_by_GO_for_" + str(name) +".fasta"]
input = [
    "pipeline/output/filtered_by_GO_cleavable_mts_for_human.fasta",
    "pipeline/output/filtered_by_GO_no_cleavable_mts_for_human.fasta",
    "pipeline/cache/filtered_proteins_by_GO_for_human.fasta"
]

#wanted_result either "absolute" or "hgt"
wanted_result = "absolute"

all_arrays = [0, 0, 0]
all_counted_instances = [0, 0, 0]
amount_of_proteins = [0, 0, 0]
for file in range(len(input)):
    panda_df = fasta_to_dataframe(input[file])
    proteomes = list(panda_df['sequence'])
    #print(proteomes)

    proteome_array = np.zeros((len(proteomes), 20), dtype=object)
    for i, proteome in enumerate(proteomes):
        proteome = proteome[1:20] 
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

    percentage = np.zeros((len(amino_acid)), dtype=object)
if wanted_result == "absolute":
    for file in range(2):
        counted_instances = all_counted_instances[file]
        cols = len(amino_acid)
        rows = len(position)
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
        cols = len(amino_acid)
        rows = len(position)
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
                    if i ==10 and j == 1:
                        print(p_value)
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
 #       print(visual_array)
 #       print(visual_array.shape)


fig, ax = plt.subplots(ncols= 2)
fig.subplots_adjust(bottom=0.5)

ax[0].set_xticks(range(len(position)), labels=position, rotation=0, rotation_mode="anchor")
ax[0].set_yticks(range(len(amino_acid)), labels=amino_acid, rotation=0, rotation_mode="anchor")
ax[0].set_yticklabels(amino_acid)
ax[0].set_xticklabels(position)

# Normalize the colormap across both graphs
vmin = min(np.min(all_arrays[0]), np.min(all_arrays[1]))
vmax = max(np.max(all_arrays[0]), np.max(all_arrays[1]))

if wanted_result == "absolute":
    cmap = "Blues" #Blues or coolwarm

if wanted_result == "hgt":
    cmap = "coolwarm" #Blues or coolwarm
pcm = ax[0].imshow(all_arrays[0], cmap=cmap, vmin=vmin, vmax=vmax)
pcm = ax[1].imshow(all_arrays[1], cmap=cmap, vmin=vmin, vmax=vmax)

ax[0].set_xlabel("Position")
ax[0].set_ylabel("Amino acids")
ax[0].set_title("Mitochondrial Proteins with MTS")

ax[1].set_xticks(range(len(position)), labels=position, rotation=0, rotation_mode="anchor")
ax[1].set_yticks(range(len(amino_acid)), labels=amino_acid, rotation=0, rotation_mode="anchor")
ax[1].set_title("Mitochondrial Proteins without MTS")
ax[1].set_xlabel("Position")
ax[1].set_ylabel("Amino acids")



fig.tight_layout(pad=3.0)
fig.colorbar(pcm, ax = ax, shrink = 0.6)
plt.show()