import pandas as pd
import requests
from tqdm import tqdm
import os
import time

# Liste der GO-Terms für ER-verwandte Lokalisationen
er_gos = {
    "GO:0005783", "GO:0005790", "GO:0005791", "GO:0009511", "GO:0016529",
    "GO:1990007", "GO:0005766", "GO:0005767", "GO:0020020", "GO:0032010",
    "GO:0036019", "GO:0042582", "GO:0042629", "GO:0043246", "GO:0044194",
    "GO:0044754", "GO:0005764", "GO:0150051", "GO:0005794", "GO:0001533",
    "GO:0042383", "GO:0097524", "GO:0120001", "GO:0005886"
}

mito_go = "GO:0005739"


def classify_location(go_terms):
    is_mito = mito_go in go_terms
    is_er = bool(er_gos & go_terms)
    if is_mito and is_er:
        return "Multiple"
    elif is_mito:
        return "GO:0005739"
    elif is_er:
        return "GO:0005783"
    else:
        return "cyto_nuclear"

def batch_query_uniprot(uniprot_ids, batch_size=500):
    """
    Holt Subcellular GO-Terms für eine Liste von UniProt-IDs in Batches.
    Gibt ein Dict mit ID → Set[GO-Terms] zurück.
    """
    base_url = "https://rest.uniprot.org/uniprotkb/search"
    headers = {"Accept": "application/json"}
    result = {}

    for i in tqdm(range(0, len(uniprot_ids), batch_size), desc="Batches"):
        batch = uniprot_ids[i:i+batch_size]
        query = " OR ".join([f"(accession:{uid})" for uid in batch])
        params = {
            "query": query,
            "fields": "accession,go",
            "format": "json",
            "size": batch_size
        }

        try:
            response = requests.get(base_url, params=params, headers=headers)
            if response.status_code != 200:
                print(f"Fehler bei Batch {i}: {response.status_code}")
                continue
            data = response.json()
            for entry in data.get("results", []):
                acc = entry["primaryAccession"]
                go_terms = {
                    ref["id"]
                    for ref in entry.get("uniProtKBCrossReferences", [])
                    if ref.get("database") == "GO"
                }
                result[acc] = go_terms
        except Exception as e:
            print(f"Fehler bei Batch {i}: {str(e)}")
        time.sleep(1)  # Optionaler Delay zur Entlastung von UniProt

    return result

def annotate_file(input_file, output_file):
    df = pd.read_csv(input_file)
    if 'ProteinID' not in df.columns:
        raise ValueError("Die CSV-Datei muss eine Spalte 'ProteinID' enthalten.")

    uniprot_ids = df['ProteinID'].dropna().unique().tolist()
    go_mapping = batch_query_uniprot(uniprot_ids)

    # Klassifizieren
    df['GO_Location'] = df['ProteinID'].apply(
        lambda pid: classify_location(go_mapping.get(pid, set()))
    )

    df.to_csv(output_file, index=False)
    print(f"✅ Datei gespeichert unter: {output_file}")

if __name__ == "__main__":
    general_working_dir = "pipeline/output/output_20250603_145910_ml_all_organisms"
    organism_names = [
    "Homo_sapiens","Mus_musculus", "Dario_rerio", "Daphnia_magna", 
    "Caenorhabditis_elegans", "Drosophila_Melanogaster", "Arabidopsis_thaliana", 
    "Physcomitrium_patens", "Chlamydomonas_reinhardtii", 
    "Candida_glabrata", "Saccharomyces_cerevisiae", "Zygosaccharomyces_rouxii"]
    organism_names = ["Arabidopsis_thaliana"]
    for name in tqdm(organism_names, leave=False):
        input_file = os.path.join(general_working_dir, name, "feature_matrix_with_go_terms.csv")
        output_file = os.path.join(general_working_dir, name, "feature_matrix_with_go_terms_alt.csv")
        
        # Load the input feature matrix
        feature_matrix1 = pd.read_csv(input_file)
        feature_matrix = feature_matrix1.copy()
        
        # Ensure 'ProteinID' column exists
        if 'protein_id' not in feature_matrix.columns:
            raise ValueError(f"Die Datei {input_file} muss eine Spalte 'protein_id' enthalten.")
        
        # Annotate GO_Term column
        uniprot_ids = feature_matrix['protein_id'].dropna().unique().tolist()
        go_mapping = batch_query_uniprot(uniprot_ids)
        
        feature_matrix['GO_Term'] = feature_matrix['protein_id'].apply(
            lambda pid: classify_location(go_mapping.get(pid, set()))
        )
        
        # Save the updated feature matrix
        feature_matrix.to_csv(output_file, index=False)
        
        # Compare the two feature matrices
        matching_rows = feature_matrix1['GO_Term'] == feature_matrix['GO_Term']
        percentage_match = (matching_rows.sum() / len(feature_matrix1)) * 100
        print(f"Percentage of matching GO_Term values for {name}: {percentage_match:.2f}%")
