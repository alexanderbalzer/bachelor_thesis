import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio import SeqIO
from scipy.stats import hypergeom
import os
from utils import transform_labels_to_names
from matplotlib.colors import TwoSlopeNorm

def fasta_to_dataframe(fasta_file):
    data = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        data.append([str(record.seq)])
    return pd.DataFrame(data, columns=['sequence'])




def run(organism_names, input_dir, cache_dir, output_dir, create_heatmap, heatmap_type, create_phylogenetic_tree, phylo_tree_type, reference):

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
        if reference == "subset":
            input = [
                os.path.join(cache_dir, f"{organism}_filtered_by_go_and_mts.fasta"),
                os.path.join(cache_dir, f"filtered_proteins_by_GO_for_{organism}.fasta")
                ]
        elif reference == "proteome":
            input = [
                os.path.join(cache_dir, f"{organism}_filtered_by_go_and_mts.fasta"),
                os.path.join(input_dir, f"{organism}.fasta")
                ]
        y = 0
        for file in input:
            data = []
            for line in SeqIO.parse(file, "fasta"):
                data.append(list(str(line.seq)))
            data = list(data)
            # limit the length of the sequences to 20 and cut the first
            for i in range(len(data)):
                data[i] = data[i][1:20]
            df = pd.DataFrame(data, columns=[list(range(19))])
            # count the instances of each amino acid per column
            df_counts = df.apply(pd.Series.value_counts).fillna(0)
            # Rearrange rows of df_counts to align with the order of amino_acid
            df_counts = df_counts.reindex(amino_acid, fill_value=0)
            # only use the first column of the dataframe
            df_counts = df_counts.iloc[:, 0]
            # transpose the dataframe
            df_counts = df_counts.transpose()
            # Convert the DataFrame to a dictionary
            count_dict = df_counts.to_dict()   
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
    ax.set_yticklabels(transform_labels_to_names(organism_names), fontstyle="italic")
    cmap = 'coolwarm'
    pcm = ax.imshow(visual_array, cmap=cmap)
    cbar = plt.colorbar(pcm, ax=ax, shrink=0.3, aspect=10, pad=0.01)
    max = np.max(visual_array)
    max = np.round(max, decimals=0)
    norm = TwoSlopeNorm(vmin=-max, vcenter=0, vmax=max)
    pcm.set_norm(norm)
    cbar.set_ticks([-max, 0, max])
    cbar.ax.set_title('HGT', pad=10)

    # Adjust the layout to prevent labels from being cut off
    plt.subplots_adjust(left=0.3)  # Increase the left margin
    plt.tight_layout()

    if create_heatmap:
        plt.savefig(os.path.join(output_dir, "heatmap.png"), dpi=300)
    # save the heatmap
    return

if __name__ == "__main__":
    # Example usage
    organism_names = ["Arabidopsis_thaliana",
    "Caenorhabditis_elegans",
    "Candida_glabrata",
    "Clavispora_lusitaniae",
    "Debaryomyces_hansenii",
    "Drosophila_Melanogaster",
    "Geotrichum_candidum",
    "human",
    "human_with_isoforms",
    "Lachancea_thermotolerans",
    "Mus_musculus",
    "Physcomitrium_patens",
    "Saccharomyces_cerevisiae",
    "Scheffersomyces_stipitis",
    "Schizosaccharomyces_pombe",
    "Yarrowia_lipolytica",
    "Zygosaccharomyces_rouxii"]
    input_dir = "pipeline/input"
    cache_dir = "pipeline/cache/cache_20250414_224234/"
    output_dir = "pipeline/output/output_20250414_224234"
    create_heatmap = True
    heatmap_type = "hgt"
    create_phylogenetic_tree = True
    phylo_tree_type = "hgt"
    reference = "proteome"
    run(organism_names, input_dir, cache_dir, output_dir, create_heatmap, heatmap_type, create_phylogenetic_tree, phylo_tree_type, reference)
