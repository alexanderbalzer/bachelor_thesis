import pandas as pd
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import os
from Bio import SeqIO
import hydrophobic_moment
from datetime import datetime
from tqdm import tqdm
import warnings
warnings.filterwarnings("ignore", category=UserWarning, module='Bio.SeqUtils.ProtParam')


def feature_matrix_per_protein(protein_sequence):
    '''
    creates a simple Feature Matrix for the given Protein
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
    get information about the MTS sequence from the MitoFates output file
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
    '''
    get the cleavage position and propbability from the signalP output file
    '''
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
Eisenberg_scale = {
    'A':  0.62,  # Alanine
    'R': -2.53,  # Arginine
    'N': -0.78,  # Asparagine
    'D': -0.90,  # Aspartic acid
    'C':  0.29,  # Cysteine
    'Q': -0.85,  # Glutamine
    'E': -0.74,  # Glutamic acid
    'G':  0.48,  # Glycine
    'H': -0.40,  # Histidine
    'I':  1.38,  # Isoleucine
    'L':  1.06,  # Leucine
    'K': -1.50,  # Lysine
    'M':  0.64,  # Methionine
    'F':  1.19,  # Phenylalanine
    'P':  0.12,  # Proline
    'S': -0.18,  # Serine
    'T': -0.05,  # Threonine
    'W':  0.81,  # Tryptophan
    'Y':  0.26,  # Tyrosine
    'V':  1.08   # Valine
}

aa_charge = {'E': -1, 'D': -1, 'K': 1, 'R': 1}




def run(organism_names, input_dir, working_dir):
    '''
    creates a Feature Matrix for the given organism
    '''
    # create a list of all the proteins in the input directory
    protein_list = []
    go_ids = ["GO:0005739", "GO:0005783"]


    for organism in tqdm(organism_names, desc="Processing organisms", position=0, leave=False):
        invalid_count = 0
        working_dir_per_organism = working_dir + "/" + organism 
        fasta = input_dir + "/" + organism + ".fasta"
        fasta_file = os.path.join(fasta)
        go_child_dir = input_dir + "/go_childs_reduced/"
        # load the go_ids
        child_dict = {}
        for go_id in go_ids:
            with open(f"{go_child_dir}{go_id}.txt", "r") as file:
                for line in file:
                    child = line.strip().split("_")[0]
                    child_dict[child] = go_id

        # Parse GO annotations
        annotation_file = f"pipeline/input/{organism}.goa"
        go_annotations = parse_go_annotations(annotation_file)

        # check if the file exists
        if not os.path.exists(fasta_file):
            print(f"File not found: {fasta_file}")
            continue
        # read the fasta file and create a DataFrame
        df = fasta_to_dataframe(fasta_file)
        # create a list of all the proteins in the fasta file
        for index, row in tqdm(df.iterrows(), total=len(df), desc=f"Processing proteins for {organism}", leave=False, position=1):
            protein_sequence = row["Sequence"]
            protein_id = row["protein_id"]
            protein_id_human = row["protein_id_human"]
            # turn protein id and sequence into a DataFrame
            feature_matrix = pd.DataFrame()
            # add the protein id to the DataFrame
            feature_matrix["protein_id"] = [protein_id]
            feature_matrix["protein_id_human"] = protein_id_human
            feature_matrix["length"] = len(protein_sequence)
            if not protein_sequence.startswith("M"):
                invalid_count += 1
                continue
            # check if the protein sequence is long enough
            if len(protein_sequence) < 12:
                invalid_count += 1
                continue
            # check if the mts sequence only contains valid amino acids
            valid_amino_acids = set("ACDEFGHIKLMNPQRSTVWY")
            if not all(aa in valid_amino_acids for aa in protein_sequence):
                invalid_count += 1
                continue
            #mts_sequence = protein_sequence[start_of_alpha_helix:start_of_alpha_helix + length_of_alpha_helix]
            # One-hot encode the second amino acid
            second_amino_acid = protein_sequence[1] if len(protein_sequence) > 1 else None
            # One-hot encode the second amino acid
            amino_acids = "ACDEFGHIKLMNPQRSTVWY"
            one_hot_encoded_second_amino_acid = {f"Second_AA_{aa}": 1 if second_amino_acid == aa else 0 for aa in amino_acids}
            feature_matrix = feature_matrix.assign(**one_hot_encoded_second_amino_acid)

            # cut the protein sequence to the length of the MTS
            cut = 40
            mts_sequence = protein_sequence[:cut]
            # if the second amino acid is A, C, T, S, V, P or G, delete the first amino acid
            if second_amino_acid in ["A", "C", "T", "S", "V", "G", "P"]:
                mts_when_huntington = mts_sequence[1:]
                # if the new first amino acid is A, C, T, S, V or G, add X as first amino acid
                if second_amino_acid in ["A", "C", "T", "S", "V", "G"]:
                    mod_mts_sequence = "X" + mts_sequence[1:]
                    alternative_mts_sequence = "XML" + mts_sequence[2:]
                else: 
                    mod_mts_sequence = mts_sequence[1:]
                    alternative_mts_sequence = mts_sequence[1:]
            # if the second amino acid is D, E, N or Q, L, I, F, Y, add X as first amino acid
            elif second_amino_acid in ["D", "E", "N", "Q", "L", "I", "F", "Y"]:
                mod_mts_sequence = "X" + mts_sequence
                alternative_mts_sequence = "XA" + mts_sequence[2:]
                mts_when_huntington = mod_mts_sequence
            else:
                mod_mts_sequence = mts_sequence
                alternative_mts_sequence = mts_sequence
                mts_when_huntington = mts_sequence
            # get the combined hydrophobicity of the first 2 amino acids of the neo n terminus
            if len(mod_mts_sequence) < 2:
                invalid_count += 1
                continue
            hyd_second_aa = Eisenberg_scale.get(mts_sequence[1], 0)
            # add the hydrophobicity of the second, third and fourth amino acid to the DataFrame
            feature_matrix["Hydrophobicity of the second amino acid"] = hyd_second_aa
            # add the hydrophobicity of the third amino acid to the DataFrame
            hyd_third_aa = Eisenberg_scale.get(mts_sequence[2], 0)
            feature_matrix["Hydrophobicity of the third amino acid"] = hyd_third_aa
            # add the hydrophobicity of the fourth amino acid to the DataFrame
            hyd_fourth_aa = Eisenberg_scale.get(mts_sequence[3], 0)
            feature_matrix["Hydrophobicity of the fourth amino acid"] = hyd_fourth_aa
            # add the helix score of the second second, third and fourth amino acid to the DataFrame
            helix_second_aa = helix_propensity.get(mts_sequence[1], 0)
            feature_matrix["Helix score of the second amino acid"] = helix_second_aa
            helix_third_aa = helix_propensity.get(mts_sequence[2], 0)
            feature_matrix["Helix score of the third amino acid"] = helix_third_aa
            helix_fourth_aa = helix_propensity.get(mts_sequence[3], 0)
            feature_matrix["Helix score of the fourth amino acid"] = helix_fourth_aa

            if mts_sequence[1] in ["A", "C", "T", "S", "V", "G", "P"] and mts_sequence[2] in ["A"]:
                feature_matrix["iMet cleavage and third amino acid is A"] = 1
            else:
                feature_matrix["iMet cleavage and third amino acid is A"] = 0



            # calculate the hydrophobic moment of the mts sequence
            hydrophobic_moment_value_mod_mts, start_best_window, length_best_window, electrostatic_help, discrimination_factor, helix_score, charge = hydrophobic_moment.run(mod_mts_sequence, verbose=False)
            # add the hydrophobic moment to the DataFrame
            feature_matrix["Hydrophobic Moment"] = hydrophobic_moment_value_mod_mts
            feature_matrix["Charge"] = charge
            feature_matrix["start_of_alpha_helix"] = start_best_window
            feature_matrix["length_of_alpha_helix"] = length_best_window
            feature_matrix["Electrostatic Help"] = electrostatic_help
            feature_matrix["Discrimination Factor"] = discrimination_factor
            '''feature_matrix["electrostatic help diff if diff nat"] = electrostatic_help - electrostatic_help_diff_nat
            feature_matrix["electrostatic help diff if huntington"] = electrostatic_help - electrostatic_help_huntington'''
            feature_matrix['helix_score'] = helix_score
            # add the protein sequence to the DataFrame
            feature_matrix["Sequence"] = [mts_sequence]
            # create a feature matrix for the protein
            feature_matrix2 = feature_matrix_per_protein(mts_sequence)
            # add the feature matrix to the feature matrix
            feature_matrix = pd.concat([feature_matrix, feature_matrix2], axis=1)
            feature_matrix['Molecular Weight'] = feature_matrix['Molecular Weight']/ cut

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
            if terms:
                # if the protein has GO terms, set the GO term to the first one
                feature_matrix['GO_Term'] = terms[0]
                # if the protein has multiple GO terms, set the GO term to "Multiple"
                if len(terms) > 1:
                    feature_matrix['GO_Term'] = "Multiple"
            #else:
            #    feature_matrix['GO_Term'] = ""
            #    terms = "Multiple"
            #    feature_matrix['GO_Term'] = terms
            #    protein_list.append(feature_matrix)
            else:
                # if the protein has no GO terms, set the GO term to cyto_nuclear
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
    # Get the absolute path to the folder this script is in
    script_dir = os.path.dirname(os.path.abspath(__file__))

    # Build the relative path to the data file
    pipeline_output_path = os.path.join(script_dir, '..', 'pipeline', 'output')
    pipeline_input_path = os.path.join(script_dir, '..', 'pipeline', 'input')

    dirs = [d for d in os.listdir(pipeline_output_path) if os.path.isdir(os.path.join(pipeline_output_path, d))]

    # Sort alphabetically and get the last one
    if dirs:
        last_dir = sorted(dirs)[-1]
        last_output_dir = os.path.join(pipeline_output_path, last_dir)
        print("Last directory:", last_output_dir)
    else:
        print("No directories found.")

    # Normalize the path (optional but good practice)
    working_dir = os.path.normpath(last_output_dir)
    input_dir = pipeline_input_path
    # Read organism names from FASTA file names
    fasta_files = [f for f in os.listdir(input_dir) if f.endswith(".fasta")]
    organism_names = [os.path.splitext(f)[0] for f in fasta_files]
    start_time = datetime.now()
    run(organism_names, input_dir, working_dir)
    end_time = datetime.now()
    run_time = end_time - start_time
    print(f"finished in: {run_time}")