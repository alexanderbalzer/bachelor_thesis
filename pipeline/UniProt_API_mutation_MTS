import requests
import json
import pandas as pd
from Bio import SeqIO
from tqdm import tqdm
import subprocess
import os


def run_perl_script(mitofates_path, input_file, flag, output_file):

    env = os.environ.copy()
    env["PATH"] = "/usr/bin:" + env["PATH"]  
    with open(output_file, "w") as output_handle:
        subprocess.run(
            ["perl", mitofates_path, input_file, flag],
            stdout=output_handle,
            stderr=subprocess.PIPE,
            env=env,
            check=True
        )




def run(organism_name, fasta_file, skip_UniProt, skip_MitoFates):
    if not skip_UniProt:
        data = []
        # Read the .fasta file using SeqIO
        with open(fasta_file, "r") as file: 
            current_entry = {}
            for line in file:
                if line.startswith(">"):
                    # Extract the protein ID from the header line
                    if current_entry:  # Save the previous entry if it exists
                        data.append(current_entry)
                    current_entry = {"Protein_ID": line.strip().split("|")[1]}
                elif line.startswith("+"):
                    current_entry["MTS_sequence"] = line.strip().split("+")[1]
                elif line.startswith("<"):
                    # Extract the sequence from the line
                    current_entry["Sequence"] = line.strip().split("<")[1]
                elif line.startswith("*"):
                    # Extract the start of MTS from the line
                    current_entry["Start_of_MTS"] = line.strip().split("*")[1]
            if current_entry:  # Append the last entry
                data.append(current_entry)

        # Create a DataFrame from the data
        df = pd.DataFrame(data)
        # Print the DataFrame
        relevant_variants = []
        for index, row in tqdm(df.iterrows()):
            protein_id = row["Protein_ID"]
            print(f"Processing protein ID: {protein_id}")
            # Extract the MTS sequence and start position
            start_of_mts = row["Start_of_MTS"]

            # URL for accessing UniProt data via API (you can choose JSON or XML format)
            url = f"https://www.ebi.ac.uk/proteins/api/variation/{protein_id}?format=json"

            # Send GET request to fetch the protein data
            response = requests.get(url)

            if response.status_code == 200:
                # Parse the response JSON
                protein_data = response.json()
                # Check if the protein entry has 'features' (mutations)
                if 'features' in protein_data:
                    variants = protein_data['features']
                    for variant in variants:
                        if int(variant.get('begin')) == int(2):  # Use .get() to avoid KeyError
                            descriptions = variant.get('association', [])  # Default to an empty list if 'descriptions' is missing
                            print("yay")
                            if descriptions:
                                disorders = []
                                # Extract the relevant information from the variant
                                for description in descriptions:
                                    disorders.append(description["name"])
                                    # Append the relevant variant information to the list
                                wildtype = variant.get('wildType')
                                mutant = variant.get('mutatedType')
                                if not mutant:
                                    continue
                                # Replace the amino acid at the position "start_of_mts" with the amino acid "mutant"
                                mutated_sequence = list(row["Sequence"])
                                mutated_sequence[int(start_of_mts) - 1] = mutant  # Adjust for 0-based indexing
                                mutated_sequence = ''.join(mutated_sequence)

                                relevant_variants.append({
                                    'Protein_ID': protein_id,
                                    'MTS_sequence': row["MTS_sequence"],
                                    'Start_of_MTS': start_of_mts,
                                    'Mutation': f"{wildtype}-{mutant}",
                                    'Identification': f"{protein_id},{wildtype}-{mutant}",
                                    'Variant_Description': disorders,
                                    'original_sequence': row["Sequence"],
                                    'Mutated_sequence': mutated_sequence,
                                    'description': variant.get('descriptions', ''),
                                    'genomic_location': variant.get('genomicLocation', ''),
                                })
            else:
                print(f"Failed to retrieve data. HTTP Status Code: {response.status_code}")
        # Convert the list of relevant variants to a DataFrame
        df_relevant_variants = pd.DataFrame(relevant_variants)
        # Save the DataFrame to a CSV file
        output_file = f"pipeline/output/output_20250507_095928/{organism_name}/MTS_variants_cache.csv"
        df_relevant_variants.to_csv(output_file, index=False)


        # save the relevant variants to a fasta file
        output_fasta_file = f"pipeline/output/output_20250507_095928/{organism_name}/{organism_name}_pathogen_variants.fasta"
        with open(output_fasta_file, "w") as fasta_file:
            for variant in relevant_variants:
                fasta_file.write(f">{variant['Protein_ID']},original\n")
                fasta_file.write(f"{variant['original_sequence']}\n")
                fasta_file.write(f">{variant['Protein_ID']},{variant['Mutation']}\n")
                fasta_file.write(f"{variant['Mutated_sequence']}\n")
    if not skip_MitoFates:
        if skip_UniProt:
            output_fasta_file = f"pipeline/output/output_20250507_095928/{organism_name}/{organism_name}_pathogen_variants.fasta"
            df_relevant_variants = pd.read_csv(f"pipeline/output/output_20250507_095928/{organism_name}/MTS_variants_cache.csv")
            relevant_variants = df_relevant_variants.to_dict(orient="records")
        # Run MitoFates as a subprocess to calculate the probability of MTS for the fasta_file
        mitofates_output_file = f"pipeline/output/output_20250507_095928/{organism_name}/MTS_probabilities.txt"
        mitofates_path = "/home/abalzer/MitoFates/MitoFates.pl"
        flag = "metazoa"
        run_perl_script(mitofates_path, output_fasta_file, flag, mitofates_output_file)
        # Read the MitoFates output file
        probability_of_mts = []
        with open(mitofates_output_file, "r") as file:
            file_without_header = file.readlines()[1:]  
            for line in file_without_header:
                fields = line.strip().split("\t")
                if len(fields) < 3:
                    continue
                protein_id = fields[0]
                probability_of_mts.append((protein_id, float(fields[1])))  # Append as a tuple

        MTS_probability_original_proteins = {}
        for identification, probability in probability_of_mts:
            protein_id, mutation_type = identification.split(",")  # Split into two variables
            if mutation_type == "original":
                MTS_probability_original_proteins[protein_id] = probability
        # Create a dictionary from the list of tuples
        probability_of_mts = dict(probability_of_mts)

        for variant in relevant_variants:
            identification = variant["Identification"]
            protein_id, mutation_type = identification.split(",")  # Split into two variables
            if mutation_type != "original":
                # Get the original protein ID from the identification
                original_protein_id = protein_id.split(",")[0]
                # Get the MTS probability for the original protein
                original_probability = MTS_probability_original_proteins.get(original_protein_id, 0.0)
                # Calculate the change in probability
                change_in_probability = f"{original_probability}->{probability_of_mts.get(identification, 0.0)}"
                delta_probability = original_probability - probability_of_mts.get(identification, 0.0)
                # Add the change in probability to the variant
                variant["Change_in_Probability_of_MTS"] = change_in_probability
                variant["Delta_Probability_of_MTS"] = delta_probability
                variant["relevant_change_in_probability"] = "*" if delta_probability > 0.1 else ""

        # Convert the list of relevant variants to a DataFrame
        df_variants = pd.DataFrame(relevant_variants)
        # Save the DataFrame to a CSV file
        output_file = f"pipeline/output/output_20250507_095928/{organism_name}/MTS_variants.csv"
        # Move the specified columns to the end of the DataFrame
        columns_to_move = ["original_sequence", "Mutated_sequence", "Identification"]
        df_variants = df_variants[[col for col in df_variants.columns if col not in columns_to_move] + columns_to_move]
        df_variants.drop(columns_to_move, axis=1, inplace=True)
        for rows in df_variants.iterrows():
            df_variants["Mutation"] = df_variants["Mutation"].str.replace("-", ">")
        # Move rows with relevant_change_in_probability = "*" to the top
        df_variants = pd.concat([
            df_variants[df_variants["relevant_change_in_probability"] == "*"],
            df_variants[df_variants["relevant_change_in_probability"] != "*"]
        ])
        df_variants.to_csv(output_file, index=False)

    if skip_MitoFates:
        df_variants = pd.read_csv(f"pipeline/output/output_20250507_095928/{organism_name}/MTS_variants.csv")
    




if __name__ == "__main__":
    organism_name = "human"  # Example organism name
    fasta_file = f"pipeline/output/output_20250507_095928/{organism_name}/MTS_sequences_{organism_name}.fasta"
    skip_UniProt = False
    skip_MitoFates = False
    run(organism_name, fasta_file, skip_UniProt, skip_MitoFates)
    
