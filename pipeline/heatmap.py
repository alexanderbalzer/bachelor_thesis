import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio import SeqIO
from scipy.stats import hypergeom
import os

def fasta_to_dataframe(fasta_file):
    data = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        data.append([str(record.seq)])
    return pd.DataFrame(data, columns=['sequence'])




def run(organism_names, cache_dir, output_dir, create_heatmap, heatmap_type, create_phylogenetic_tree, phylo_tree_type):

    amino_acid = np.array(["D", "E", "N", "Q", "Y", "H", "K", "R", "M", "L", "F", "I", "W", "S", "A", "T", "C", "P", "G", "V"])

    if create_phylogenetic_tree:
        if phylo_tree_type == "absolute":
            save_subset_array_for_phylogenetic_tree = True
            save_HGT_array_for_phylogenetic_tree = False
        elif phylo_tree_type == "hgt":
            save_subset_array_for_phylogenetic_tree = False
            save_HGT_array_for_phylogenetic_tree = True
        else:
            raise ValueError("Invalid phylogenetic tree type. Choose 'absolute' or 'hgt'.")

    dictlist = []
    reference_array = np.zeros((len(organism_names), len(amino_acid)), dtype=float)
    subset_array = np.zeros((len(organism_names), len(amino_acid)), dtype=float)
    visual_array = np.zeros((len(organism_names), len(amino_acid)), dtype=float)
    for z, organism in enumerate(organism_names):
        input = [
            os.path.join(output_dir, f"{organism}_filtered_by_go_and_mts.fasta"),
            os.path.join(cache_dir, f"filtered_proteins_by_GO_for_{organism}.fasta")
            ]
        y = 0
        for file in input:
            panda_df = fasta_to_dataframe(file)
            proteome = list(panda_df['sequence'])
            amount_of_proteins = len(proteome)

            proteome_as_array = np.zeros((amount_of_proteins, 20), dtype=object)
            for i, amino_acids_in_proteome in enumerate(proteome):
                sequence = amino_acids_in_proteome[1:20]
                for j in range(19):
                    if j < len(sequence):
                        proteome_as_array[i, j] = sequence[j]
                    else:
                        proteome_as_array[i, j] = "X"
                count_dict = {}
                for row in range(amount_of_proteins):
                    value = proteome_as_array[row, 0]
                    if (value) in count_dict:
                        count_dict[value] += 1
                    else:
                        count_dict[value] = 1
            dictlist.append(count_dict)
            if y == 0:
                for i in range(len(amino_acid)):
                    subset_array[z, i] = count_dict.get(amino_acid[i], 0)
            if y == 1:
                for i in range(len(amino_acid)):
                    reference_array[z, i] = count_dict.get(amino_acid[i], 0)
            y += 1

    # Calculate HGT scores:
        for i in range(len (organism_names)):
            M = np.sum(reference_array[i])
            N = np.sum(subset_array[i])

            for j in range(len(amino_acid)):
                x = subset_array[i, j]
                n = reference_array[i, j]
                p_value = hypergeom.pmf(x, M, n, N)
                abs_log_val = abs(np.log10(p_value))
                f_obs = x / N if N != 0 else 0
                f_exp = n / M if M != 0 else 0
                result = 0
                if f_obs < f_exp:
                    result = -1 * abs_log_val
                else:
                    result = abs_log_val
                visual_array[i, j] = result
    # Save the visual array to a file
    if save_HGT_array_for_phylogenetic_tree:
        np.save(os.path.join(cache_dir, "phyl_tree_array.npy"), visual_array)

    # normalize the subset array
    for i in range(len(organism_names)):
        for j in range(len(amino_acid)):
            if subset_array[i, j] != 0:
                subset_array[i, j] = subset_array[i, j] / np.sum(subset_array[i])
            else:
                subset_array[i, j] = 0
    # Save the subset array to a file
    #print(subset_array)
    if save_subset_array_for_phylogenetic_tree:
        np.save(os.path.join(cache_dir, "phyl_tree_array.npy"), subset_array)


    if create_heatmap:
        if heatmap_type == "absolute":
            visual_array = subset_array
        elif heatmap_type == "hgt":
            visual_array = visual_array
        else:
            raise ValueError("Invalid heatmap type. Choose 'absolute' or 'hgt'.")
    fig, ax = plt.subplots()
    ax.set_xticks(np.arange(len(amino_acid)), amino_acid)
    ax.set_yticks(np.arange(len(organism_names)), organism_names)
    ax.set_xticklabels(amino_acid)
    ax.set_yticklabels(organism_names)
    cmap = 'coolwarm'
    pcm = ax.imshow(visual_array, cmap=cmap)
    plt.colorbar(pcm, ax=ax, shrink=0.8)
    plt.title("HGT scores")
    plt.show()
    # save the heatmap
    if create_heatmap:
        plt.savefig(os.path.join(output_dir, "heatmap.png"), dpi=300)
        plt.close()
    return

