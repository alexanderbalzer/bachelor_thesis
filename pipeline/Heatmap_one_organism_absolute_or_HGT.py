import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio import SeqIO
from scipy.stats import hypergeom
from matplotlib.colors import TwoSlopeNorm
import os

def fasta_to_dataframe(fasta_file):
    data = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        data.append([str(record.seq)])
    return pd.DataFrame(data, columns=['sequence'])


def count_instances_at_positions(array):
    rows, cols = array.shape
    dictlist = []
    for col in range(0,20): 
        count_dict = {}
        for row in range(rows):
            value = array[row, col]
            if (value) in count_dict:
                count_dict[value] += 1
            else:
                count_dict[value] = 1
        dictlist.append(count_dict)
    return dictlist



def run(organism_names, input_dir, output_dir, heatmap_type):
    """
    Main function to run the pipeline. Takes a list of organisms, loads their fasta files, counts instances of amino acids at each position,
    and generates heatmaps based on the specified heatmap type (absolute or HGT).
    Args:
        organism_names (list): List of organism names.
        input_dir (str): Directory containing the input FASTA files.
        output_dir (str): Directory to save the output heatmaps.
        heatmap_type (str): Type of heatmap to generate ("absolute" or "hgt").
    Outputs:
        Generates heatmaps for each organism and saves them in the specified output directory.
    """
    for i, name in enumerate(organism_names, start=1):
        print(f"Processing organism: {name}")
        output_dir_per_organism = output_dir + "/" + name 
        position = np.array(["2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20"])
        amino_acid = np.array(["D", "E", "N", "Q", "Y", "H", "K", "R", "M", "L", "F", "I", "W", "S", "A", "T", "C", "P", "G", "V"])

        # Specify the input files for cleavable, non-cleavable MTS and reference
        input = [
            output_dir_per_organism + "/" + name + "_filtered_by_GO_cleavable_mts.fasta",
            output_dir_per_organism + "/" + name + "_filtered_proteins_by_GO_noncleavable_mts.fasta",
            os.path.join(input_dir, f"{name}.fasta")
            ]
        # Check if the input files exist
        for file in input:
            if not os.path.exists(file):
                print(f"Input file not found: {file}")
                continue
        # setup the arrays
        all_arrays = [0, 0]
        all_counted_instances = [0, 0, 0]
        amount_of_proteins = [0, 0, 0]
        subset_protein_count = {}
        # run the three input files
        for file in range(len(input)):
            panda_df = fasta_to_dataframe(input[file])
            subset_protein_count[file] = len(panda_df)
            proteins = list(panda_df['sequence'])

            # Create a 2D array to store the first 20 amino acids of each protein
            proteome_array = np.zeros((len(proteins), 20), dtype=object)
            for i, protein in enumerate(proteins):
                protein = protein[1:21] 
                for j in range(20):
                    if j < len(protein):
                        proteome_array[i, j] = protein[j]
                    else:
                        proteome_array[i, j] = 'X'

            # Count the instances of each amino acid at each position and store in a dictionary
            counted_instances = count_instances_at_positions(proteome_array)
            rows, cols = proteome_array.shape
            amount_of_proteins[file] = rows
            all_counted_instances[file] = counted_instances

        # Count the instances of each amino acid at each position in the two subsets
        if heatmap_type == "absolute":
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
                
        # Calculate the HGT score for each amino acid in the two subsets 
        if heatmap_type == "hgt":
            whole_set = all_counted_instances[2]
            for file in range(2):
                counted_instances = all_counted_instances[file]
                cols = len(amino_acid)
                rows = len(position)
                visual_array = np.zeros((cols, rows))
                for i in range(rows):
                    for j in range(cols):
                        # Calculate the HGT score for each amino acid at each position using the hypergeometric distribution                        
                        if amino_acid[j] in counted_instances[i]:
                            x = counted_instances[i][amino_acid[j]] #amount of amino acid in the subset
                            n = whole_set[i][amino_acid[j]] #amount of amino acid in the whole set
                            M = amount_of_proteins[2] #amount of all amino acids in the whole set
                            N = amount_of_proteins[file] #amount of all amino acids in the subset
                            p_value = hypergeom.sf(x, M, n, N)
                            p_value_over = hypergeom.sf(x, M, n, N) # p-value for the upper tail
                            p_value_under = hypergeom.cdf(x, M, n, N) # p-value for the lower tail
                            p_value = min(p_value_over, p_value_under) # take the lower value
                            abs_log_val = abs(np.log10(p_value))
                            f_obs = x / N
                            f_exp = n / M
                            result = 0
                            if f_obs < f_exp:
                                result = -1 * abs_log_val
                            else:
                                result = abs_log_val
                            visual_array[j, i] = result
                all_arrays[file]= visual_array


        # Create the heatmap
        # Set the color map based on the heatmap type
        if heatmap_type == "absolute":
            cmap = "Blues" 
        if heatmap_type == "hgt":
            cmap = "RdBu_r" 
        fig, ax = plt.subplots(ncols= 2)
        fig.subplots_adjust(bottom=0.5)
        ax[0].set_xticks(range(len(position)), labels=position, rotation=0, rotation_mode="anchor", fontsize=7)
        ax[0].set_yticks(range(len(amino_acid)), labels=amino_acid, rotation=0, rotation_mode="anchor", fontsize=7)
        ax[0].set_yticklabels(amino_acid)
        ax[0].set_xticklabels(position)
        pcm1 = ax[0].imshow(all_arrays[0], cmap=cmap)
        pcm2 = ax[1].imshow(all_arrays[1], cmap=cmap)
        ax[0].set_ylabel("")
        ax[0].set_title(f"Mitochondrial proteins with MTS (n= {subset_protein_count[0]})", fontsize=7)
        ax[1].set_xticks(range(len(position)), labels=position, rotation=0, rotation_mode="anchor", fontsize=7)
        ax[1].set_yticks(range(len(amino_acid)), labels=amino_acid, rotation=0, rotation_mode="anchor", fontsize=7)
        ax[1].set_title(f"Mitochondrial proteins without MTS (n= {subset_protein_count[1]})", fontsize=7)
        ax[1].set_ylabel(" ")
        ax[0].set_xlabel("Position in N-terminal sequences", fontsize=7)
        ax[1].set_xlabel("Position in N-terminal sequences", fontsize=7)
        fig.tight_layout(pad=3.0)
        max = 20
        norm = TwoSlopeNorm(vmin=-max, vcenter=0, vmax=max)
        pcm1.set_norm(norm)
        pcm2.set_norm(norm)
        cbar = plt.colorbar(pcm1, ax=ax, shrink=0.3, aspect=10, pad=0.01)
        cbar.set_ticks(np.linspace(-max, max, num=3))
        cbar.set_ticklabels([str(-max), "0", str(max)], fontsize=7)
        cbar.ax.set_title('HGT', pad=5, fontsize=7)
        plt.savefig(os.path.join(output_dir_per_organism, f"heatmap_{name}_{heatmap_type}.png"), dpi=300)
        plt.close(fig)

if __name__ == "__main__":
    # Define the names of the organisms
    organism_names = [
    "Homo_sapiens","Mus_musculus", "Rattus_norvegicus", "Dario_rerio",
    "Caenorhabditis_elegans", "Drosophila_Melanogaster", "Arabidopsis_thaliana", 
    "Saccharomyces_cerevisiae"]
    
    output_dir = "pipeline/output/output_20250616_161709"
    input_dir = "pipeline/input"
    heatmap_type = "hgt"  # or "hgt"

    run(organism_names, input_dir, output_dir, heatmap_type)
    