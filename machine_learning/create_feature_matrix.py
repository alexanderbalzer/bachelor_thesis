import pandas as pd
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import os
from Bio import SeqIO
import hydrophobic_moment
from datetime import datetime
from tqdm import tqdm
import add_go_terms_to_df
from goatools.obo_parser import GODag
from goatools.base import download_go_basic_obo
import numpy as np
import re
import warnings
from sklearn.preprocessing import OneHotEncoder
warnings.filterwarnings("ignore", category=UserWarning, module='Bio.SeqUtils.ProtParam')


def feature_matrix_per_protein(protein_sequence):
    '''
    creates a Feature Matrix for the given Protein
    '''
    X = ProteinAnalysis(protein_sequence)
    # create a dictionary with the features
    features = {
        'molecular_weight': X.molecular_weight(),
        'L and A percentage': X.amino_acids_percent['L'] if 'L' in X.amino_acids_percent else 0 + X.amino_acids_percent['A'] if 'A' in X.amino_acids_percent else 0,
        'R percentage': X.amino_acids_percent['R'] if 'R' in X.amino_acids_percent else 0,
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
    flat['Molecular Weight'] = features['molecular_weight']
    flat['Leucine_and_Alanine_percentage'] = features['L and A percentage']
    flat['Arginine_percentage'] = features['R percentage']
    flat['Aromaticity'] = features['Aromaticity']
    flat['Instability Index'] = features['Instability Index']
    flat['Instability Index'] = features['Instability Index']
    flat['Isoelectric Point'] = features['Isoelectric Point']
    flat['Molecular Weight'] = features['Molecular Weight']
    flat['SecStr_Helix'] = features['Secondary Structure Fraction'][0]
    flat['SecStr_Sheet'] = features['Secondary Structure Fraction'][1]
    flat['Hydrophobicity'] = features['Hydrophobicity']

    # Store as DataFrame
    df = pd.DataFrame([flat])
    return df


def fasta_to_dataframe(fasta_file):
    sequences = []
    protein_id = []
    protein_id_human = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequences.append(str(record.seq))
        protein_id_uniprot = record.id.split("|")[1]
        protein_id_human.append(record.id.split("|")[2])
        protein_id.append(str(protein_id_uniprot))
    
    df = pd.DataFrame({
        "Sequence": sequences,
        "protein_id": protein_id,
        "protein_id_human": protein_id_human
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
            tom20_motive = fields[7]
            if tom20_motive == "-":
                tom20_motive = 0
            else:
                tom20_motive = 1
            start_of_alpha_helix = position_of_MTS.strip().split("-")[0]
            length_of_alpha_helix = float(position_of_MTS.strip().split("-")[1]) - float(position_of_MTS.strip().split("-")[0])
            length_of_MTS = fields[3].split("(")[0]
            if length_of_MTS.isnumeric():
                length_of_MTS = float(length_of_MTS)
            else:
                length_of_MTS = None
            mitofates_cleavage_probability = fields[1]
            return float(start_of_alpha_helix), float(length_of_alpha_helix), length_of_MTS, float(mitofates_cleavage_probability), float(tom20_motive)
        
def get_signalP_infos(working_dir, name, wanted_id):
    with open(working_dir + "/"+ name + "/processed_signalp_results.txt", "r") as file:
        lines = file.readlines()
        for i, line in enumerate(lines):
            if line.startswith("#") or i == 0:
                continue
            fields = line.strip().split("\t")
            protein_id = fields[0]  # protein ID (column 1 in SignalP)
            if protein_id != wanted_id:
                continue
            if len(fields) > 4:
                predicted_cleavage_position = fields[4].strip().split(":")[1].split("-")[0].strip()
            else:
                predicted_cleavage_position = None
            signalP_cleavage_probability = fields[3]
            if predicted_cleavage_position == None:
                return None, float(signalP_cleavage_probability)
            else:
                return float(predicted_cleavage_position), float(signalP_cleavage_probability)

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


helix_propensity = {
    'A': 1.45, 'C': 0.77, 'D': 0.98, 'E': 1.53,
    'F': 1.12, 'G': 0.53, 'H': 1.24, 'I': 1.00,
    'K': 1.07, 'L': 1.34, 'M': 1.20, 'N': 0.73,
    'P': 0.59, 'Q': 1.17, 'R': 0.79, 'S': 0.79,
    'T': 0.82, 'V': 1.14, 'W': 1.14, 'Y': 0.61
}

def helix_score(seq):
    if len(seq) < 10:
        return 0  # Return 0 if the sequence is shorter than the frame length
    max_score = 0
    for i in range(len(seq) - 9):  # Iterate over all possible frames of length 10
        frame = seq[i:i + 10]
        frame_score = sum(helix_propensity.get(res.upper(), 0) for res in frame) / len(frame)
        max_score = max(max_score, frame_score)
    return max_score

def classify_natc_substrate(second_as):
    if second_as in ["A", "C", "T", "S", "V", "G"]:
        return "NatA/D"
    elif second_as in ["D", "E", "N", "Q"]:
        return "NatB"
    elif second_as in ["L", "I", "F", "Y", "K"]:
        return "NatC/E"
    else:
        return "Other"
    


def run(organism_names, input_dir, working_dir):
    '''
    creates a Feature Matrix for the given organism
    '''
    # create a list of all the proteins in the input directory
    protein_list = []
    go_ids2 = [
            #Go terms for: ER, Golgi, Ribosome, Mitochondria, Nucleus, Lysosome, Cell membrane, Cytoplasm
            "GO:0005783",
            "GO:0005794",
            "GO:0005840",
            "GO:0005739",
            "GO:0005634",
            "GO:0005764",
            "GO:0005886",
            "GO:0005829"
            ]
    go_ids = ["GO:0005739", "GO:0005783"]

    for organism in organism_names:
        invalid_count = 0
        working_dir_per_organism = working_dir + "/" + organism 
        fasta = working_dir_per_organism + "/" + "filtered_proteins_by_GO_for_" + organism + ".fasta"
        fasta_file = os.path.join(fasta)
        go_child_dir = working_dir_per_organism + "/go_childs_reduced/"
        # load the go_ids
        child_dict = {}
        for go_id in go_ids:
            with open(f"{go_child_dir}{go_id}.txt", "r") as file:
                for line in file:
                    child = line.strip().split("_")[0]
                    child_dict[child] = go_id

        # Parse GO annotations
        annotation_file = "pipeline/input/Homo_sapiens.goa"
        go_annotations = parse_go_annotations(annotation_file)

        # check if the file exists
        if not os.path.exists(fasta_file):
            print(f"File not found: {fasta_file}")
            continue
        # read the fasta file and create a DataFrame
        df = fasta_to_dataframe(fasta_file)
        # create a list of all the proteins in the fasta file
        for index, row in tqdm(df.iterrows(), total=len(df), desc=f"Processing proteins for {organism}"):
            protein_sequence = row["Sequence"]
            protein_id = row["protein_id"]
            protein_id_human = row["protein_id_human"]
            # turn protein id and sequence into a DataFrame
            feature_matrix = pd.DataFrame()
            # add the protein id to the DataFrame
            feature_matrix["protein_id"] = [protein_id]
            feature_matrix["protein_id_human"] = protein_id_human
            protein_id = row["protein_id"]
            start_of_alpha_helix, length_of_alpha_helix, mpp_cleavage_pos, mitofates_cleavage_probability, tom20_motive = get_mitoFates_infos(working_dir, organism, protein_id)
            feature_matrix["start_of_alpha_helix"] = start_of_alpha_helix
            feature_matrix["length_of_alpha_helix"] = length_of_alpha_helix
            feature_matrix["MPP_cleavage_position"] = mpp_cleavage_pos
            feature_matrix["mitofates_cleavage_probability"] = mitofates_cleavage_probability
            feature_matrix["tom20_motive"] = tom20_motive
            # get the signalP information
            spi_cleavage_pos, signalP_cleavage_probability = get_signalP_infos(working_dir, organism, protein_id)
            if spi_cleavage_pos == None:
                signalP_cleavage_probability = 0
            else:
                signalP_cleavage_probability = 1
            feature_matrix["signalP_cleavage_probability"] = signalP_cleavage_probability
            # mts sequence is the sequence specified by start and length of alpha helix
            start_of_alpha_helix = int(start_of_alpha_helix) -1
            length_of_alpha_helix = int(length_of_alpha_helix)
            if not protein_sequence.startswith("M"):
                invalid_count += 1
                continue
            # check if the mts sequence only contains valid amino acids
            valid_amino_acids = set("ACDEFGHIKLMNPQRSTVWY")
            if not all(aa in valid_amino_acids for aa in protein_sequence):
                invalid_count += 1
                continue
            cleavage_pos = np.min([mpp_cleavage_pos, spi_cleavage_pos]) if spi_cleavage_pos is not None else mpp_cleavage_pos
            if len(protein_sequence) < int(cleavage_pos):
                invalid_count += 1
                continue
            #mts_sequence = protein_sequence[start_of_alpha_helix:start_of_alpha_helix + length_of_alpha_helix]
            mts_sequence = protein_sequence[0:int(cleavage_pos)]
            # calculate the hydrophobic moment of the mts sequence
            hydrophobic_moment_value = hydrophobic_moment.run(mts_sequence, verbose=False)
            # add the hydrophobic moment to the DataFrame
            feature_matrix["Hydrophobic Moment"] = hydrophobic_moment_value
            # One-hot encode the second amino acid
            second_amino_acid = protein_sequence[1] if len(protein_sequence) > 1 else None
            amino_acids = "ACDEFGHIKLMNPQRSTVWY"
            # Classify the NAT substrate
            natc_substrate = classify_natc_substrate(second_amino_acid)
            #feature_matrix["NAT_Substrate"] = natc_substrate

            # One-hot encode the NAT substrate classification
            natc_classes = ["NatA/D", "NatB", "NatC/E", "Other"]
            one_hot_encoded_natc = {f"NAT_{cls}": 1 if natc_substrate == cls else 0 for cls in natc_classes}
            feature_matrix = feature_matrix.assign(**one_hot_encoded_natc)
            # cut the protein sequence to the length of the MTS
            protein_sequence = protein_sequence[:90]
            # add the protein sequence to the DataFrame
            feature_matrix["Sequence"] = [protein_sequence]
            # create a feature matrix for the protein
            feature_matrix2 = feature_matrix_per_protein(protein_sequence)
            # add the feature matrix to the feature matrix
            feature_matrix = pd.concat([feature_matrix, feature_matrix2], axis=1)
            feature_matrix['Molecular Weight'] = feature_matrix['Molecular Weight']/ cleavage_pos
            feature_matrix['helix_score'] = helix_score(mts_sequence)

            # add the GO terms to the feature matrix

            # Add filtered GO terms rowwise to the feature matrix
            filtered_terms = []
            if protein_id in go_annotations:
                for term in go_annotations[protein_id]:
                    if term in child_dict:
                        filtered_terms.append(child_dict[term])
                # Remove duplicate GO terms in each row
                filtered_terms = list(set(filtered_terms))
                # if a protein has multiple go terms, set the filtered terms empty
                terms = filtered_terms
            else:
                feature_matrix['GO_Term'] = ""
            '''if "GO:0005739" in terms:
                if mitofates_cleavage_probability < 0.5:
                    terms = ["GO:0005739_no_cleavable_mts" if t == "GO:0005739" else t for t in terms]
                elif mitofates_cleavage_probability >= 0.5:
                    terms = ["GO:0005739_cleavable_mts" if t == "GO:0005739" else t for t in terms]'''
            if len(terms) > 1:
                '''for term in terms:
                    feature_matrix_temp = feature_matrix.copy()
                    feature_matrix_temp['GO_Term'] = term
                    protein_list.append(feature_matrix_temp)'''
                terms = "Multiple"
                feature_matrix['GO_Term'] = terms
                protein_list.append(feature_matrix)
            else:
                # if the protein has no GO terms, set the GO term to cyto_nuclear
                if len(terms) == 0:
                    terms = "cyto_nuclear"
                feature_matrix['GO_Term'] = terms
                # append the feature matrix to the list of proteins
                protein_list.append(feature_matrix)

            # append the feature matrix to the list of proteins
        all_proteins_df = pd.concat(protein_list, ignore_index=True)
        # save the DataFrame to a csv file
        output_file = os.path.join(working_dir_per_organism, "feature_matrix_with_go_terms.csv")
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