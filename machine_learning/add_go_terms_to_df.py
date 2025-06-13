import pandas as pd
import requests
from tqdm import tqdm
import os

# Liste der GO-Terms für ER-verwandte Lokalisationen
er_gos = {
    "GO:0005783", "GO:0005790", "GO:0005791", "GO:0009511", "GO:0016529",
    "GO:1990007", "GO:0005766", "GO:0005767", "GO:0020020", "GO:0032010",
    "GO:0036019", "GO:0042582", "GO:0042629", "GO:0043246", "GO:0044194",
    "GO:0044754", "GO:0005764", "GO:0150051", "GO:0005794", "GO:0001533",
    "GO:0042383", "GO:0097524", "GO:0120001", "GO:0005886"
}

mito_go = "GO:0005739"

def get_go_terms(uniprot_id):
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"
    try:
        response = requests.get(url)
        if response.status_code != 200:
            return set()
        data = response.json()
        refs = data.get("uniProtKBCrossReferences", [])
        return {ref["id"] for ref in refs if ref.get("database") == "GO"}
    except Exception:
        return set()

def classify_location(go_terms):
    is_mito = mito_go in go_terms
    is_er = bool(er_gos & go_terms)

    if is_mito and is_er:
        return "Multiple"
    elif is_mito:
        return 'GO:0005739'
    elif is_er:
        return "GO:0005783"
    else:
        return "cyto_nuclear"

def annotate_file(input_file, output_file):
    df = pd.read_csv(input_file)
    if 'ProteinID' not in df.columns:
        raise ValueError("Die CSV-Datei muss eine Spalte 'ProteinID' enthalten.")

    tqdm.pandas(desc="Annotating proteins")
    df['GO_Location'] = df['ProteinID'].progress_apply(lambda pid: classify_location(get_go_terms(pid)))

    df.to_csv(output_file, index=False)
    print(f"✅ Ergebnis gespeichert unter: {output_file}")

if __name__ == "__main__":
    general_working_dir = "pipeline/output/output_20250603_145910_ml_all_organisms"
    organism_names = [
    "Homo_sapiens","Mus_musculus", "Dario_rerio", "Daphnia_magna", 
    "Caenorhabditis_elegans", "Drosophila_Melanogaster", "Arabidopsis_thaliana", 
    "Physcomitrium_patens", "Chlamydomonas_reinhardtii", 
    "Candida_glabrata", "Saccharomyces_cerevisiae", "Zygosaccharomyces_rouxii"]
    organism_names = ["Homo_sapiens"]
    for name in organism_names:
        input_file = os.path.join(general_working_dir, name, "feature_matrix_with_go_terms")
        output_file = os.path.join(general_working_dir, name, "feature_matrix_with_go_terms_alt")
        feature_matrix1 = pd.read_csv(input_file)
        feature_matrix = feature_matrix1.copy()
        for i, row in feature_matrix.iterrows():
            protein_id = row['ProteinID']
            go_terms = get_go_terms(protein_id)
            location = classify_location(go_terms)
            feature_matrix.at[i, 'GO_Term'] = location
        feature_matrix.to_csv(output_file, index=False)
        # Compare the two feature matrices
        matching_rows = feature_matrix1['GO_Term'] == feature_matrix['GO_Term']
        percentage_match = (matching_rows.sum() / len(feature_matrix1)) * 100
        print(f"Percentage of matching GO_Term values: {percentage_match:.2f}%")
