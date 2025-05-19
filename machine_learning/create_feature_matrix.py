import pandas as pd
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import os
from Bio import SeqIO
import hydrophobic_moment
from datetime import datetime


def feature_matrix_per_protein(protein_sequence):
    '''
    creates a Feature Matrix for the given Protein
    '''
    X = ProteinAnalysis(protein_sequence)
    # create a dictionary with the features
    features = {
        'Amino Acids Percent': X.get_amino_acids_percent(),
        'molecular_weight': X.molecular_weight(),
        'Aromaticity': X.aromaticity(),
        'Instability Index': X.instability_index(),
        'Isoelectric Point': X.isoelectric_point(),
        'Molecular Weight': X.molecular_weight(),
        'Secondary Structure Fraction': X.secondary_structure_fraction(),
        'Hydrophobicity': X.gravy(),
        }
    # create a DataFrame with the features
    # Flatten
    flat = {}
    flat.update({f'AA_{k}': v for k, v in features['Amino Acids Percent'].items()})
    flat['Molecular Weight'] = features['molecular_weight']
    flat['Aromaticity'] = features['Aromaticity']
    flat['Instability Index'] = features['Instability Index']
    flat['Instability Index'] = features['Instability Index']
    flat['Isoelectric Point'] = features['Isoelectric Point']
    flat['Molecular Weight'] = features['Molecular Weight']
    flat['SecStr_Helix'] = features['Secondary Structure Fraction'][0]
    flat['SecStr_Sheet'] = features['Secondary Structure Fraction'][1]
    flat['SecStr_Turn'] = features['Secondary Structure Fraction'][2]
    flat['Hydrophobicity'] = features['Hydrophobicity']

    # Store as DataFrame
    df = pd.DataFrame([flat])
    return df


def fasta_to_dataframe(fasta_file):
    sequences = []
    protein_id = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequences.append(str(record.seq))
        protein_id_uniprot = record.id.split("|")[1]
        protein_id.append(str(protein_id_uniprot))
    
    df = pd.DataFrame({
        "Sequence": sequences,
        "protein_id": protein_id
    })
    return df

def get_mitoFates_infos(working_dir, name, wanted_id):
    '''
    get the length of the MTS sequence from the MitoFates output
    '''
    with open(working_dir + "/"+ name + "/mitofates_for_" + name + ".cgi", "r") as file:
        lines = file.readlines()
        for i, line in enumerate(lines):
            if line.startswith("!") or i == 0:  # Skip header line
                continue
            fields = line.strip().split("\t")
            protein_id = fields[0].strip().split("|")[1]  # protein ID (column 1 in MitoFates)
            if protein_id != wanted_id:
                continue
            position_of_MTS = fields[6]
            start_of_alpha_helix = position_of_MTS.strip().split("-")[0]
            length_of_alpha_helix = float(position_of_MTS.strip().split("-")[1]) - float(position_of_MTS.strip().split("-")[0])
            length_of_MTS = fields[3].split("(")[0]
            cleavage_probability = fields[1]
            return start_of_alpha_helix, length_of_alpha_helix, length_of_MTS, cleavage_probability

def parse_go_annotations(annotation_file):
    """
    Parse a GO annotation file to create a dictionary
    mapping protein IDs and their associated GO terms.
    """
    go_annotation = {}
    with open(annotation_file, "r") as file:
        for line in file:
            if line.startswith("!"):  
                continue
            fields = line.strip().split("\t")
            protein_id = fields[1]  #protein ID (column 2 in GAF)
            function = fields[8]  #function type (column 9 in GAF)
            if function != "C":  #filter for location information
                continue
            go_term = fields[4]     #GO term (column 5 in GAF)
            if protein_id not in go_annotation:
                go_annotation[protein_id] = []
            go_annotation[protein_id].append(go_term)
    return go_annotation


def run(organism_names, input_dir, working_dir):
    '''
    creates a Feature Matrix for the given organism
    '''
    # create a list of all the proteins in the input directory
    protein_list = []
    for organism in organism_names:
        invalid_count = 0
        working_dir_per_organism = working_dir + "/" + organism 
        fasta = working_dir_per_organism + "/" + "filtered_proteins_by_GO_for_" + organism + ".fasta"
        fasta_file = os.path.join(fasta)
        # check if the file exists
        if not os.path.exists(fasta_file):
            print(f"File not found: {fasta_file}")
            continue
        # read the fasta file and create a DataFrame
        df = fasta_to_dataframe(fasta_file)
        # create a list of all the proteins in the fasta file
        for index, row in df.iterrows():
            print(f"Processing protein: {index+1}/{len(df)}")
            protein_sequence = row["Sequence"]
            protein_id = row["protein_id"]
            # turn protein id and sequence into a DataFrame
            feature_matrix = pd.DataFrame()
            # add the protein id to the DataFrame
            feature_matrix["protein_id"] = [protein_id]
            protein_id = row["protein_id"]
            start_of_alpha_helix, length_of_alpha_helix, length_of_MTS, cleavage_probability = get_mitoFates_infos(working_dir, organism, protein_id)
            feature_matrix["start_of_alpha_helix"] = start_of_alpha_helix
            feature_matrix["length_of_alpha_helix"] = length_of_alpha_helix
            feature_matrix["length_of_MTS"] = length_of_MTS
            feature_matrix["cleavable_mts"] = cleavage_probability
            # mts sequence is the sequence specified by start and length of alpha helix
            start_of_alpha_helix = int(start_of_alpha_helix) -1
            length_of_alpha_helix = int(length_of_alpha_helix)
            if not protein_sequence.startswith("M"):
                invalid_count += 1
                continue
            try:
                length_of_MTS = float(length_of_MTS)
            except ValueError:
                continue
            # check if the mts sequence only contains valid amino acids
            valid_amino_acids = set("ACDEFGHIKLMNPQRSTVWY")
            if not all(aa in valid_amino_acids for aa in protein_sequence):
                print(f"Skipping protein {protein_id} due to invalid amino acids")
                invalid_count += 1
                continue
            if len(protein_sequence) < int(length_of_MTS):
                print(f"Skipping protein {protein_id} due to length of MTS")
                invalid_count += 1
                continue
            #mts_sequence = protein_sequence[start_of_alpha_helix:start_of_alpha_helix + length_of_alpha_helix]
            mts_sequence = protein_sequence[0:int(length_of_MTS)]
            # calculate the hydrophobic moment of the mts sequence
            hydrophobic_moment_value = hydrophobic_moment.run(mts_sequence)
            # add the hydrophobic moment to the DataFrame
            feature_matrix["Hydrophobic Moment"] = hydrophobic_moment_value
            # One-hot encode the second amino acid
            second_amino_acid = protein_sequence[1] if len(protein_sequence) > 1 else None
            amino_acids = "ACDEFGHIKLMNPQRSTVWY"
            one_hot_encoded = {f"Second_AA_{aa}": 1 if second_amino_acid == aa else 0 for aa in amino_acids}
            feature_matrix = feature_matrix.assign(**one_hot_encoded)
            # cut the protein sequence to the length of the MTS
            protein_sequence = protein_sequence[:int(length_of_MTS)]
            # add the protein sequence to the DataFrame
            feature_matrix["Sequence"] = [protein_sequence]
            # create a feature matrix for the protein
            feature_matrix2 = feature_matrix_per_protein(protein_sequence)
            # add the feature matrix to the feature matrix
            feature_matrix = pd.concat([feature_matrix, feature_matrix2], axis=1)

            # append the feature matrix to the list of proteins
            protein_list.append(feature_matrix)
        # create a DataFrame with all the proteins
        all_proteins_df = pd.concat(protein_list, ignore_index=True)
        # save the DataFrame to a csv file
        output_file = os.path.join(working_dir_per_organism, "feature_matrix.csv")
        all_proteins_df.to_csv(output_file, index=False)
        print(f"Feature matrix saved to {output_file}")
        print(f"Invalid proteins skipped: {invalid_count}")
    
if __name__ == "__main__":
    # Example usage
    '''organism_names = [
    "Homo_sapiens", "Homo_sapiens_isoforms", "Mus_musculus", "Dario_rerio", "Daphnia_magna", 
    "Caenorhabditis_elegans", "Drosophila_Melanogaster", "Arabidopsis_thaliana", 
    "Physcomitrium_patens", "Chlamydomonas_reinhardtii", 
    "Candida_glabrata", "Saccharomyces_cerevisiae", "Zygosaccharomyces_rouxii"]'''
    organism_names = ["Homo_sapiens"]
    working_dir = "pipeline/output/output_20250519_142700_machine_learning_human"
    input_dir = "pipeline/input"
    start_time = datetime.now()
    run(organism_names, input_dir, working_dir)
    end_time = datetime.now()
    run_time = end_time - start_time
    print(f"finished in: {run_time}")