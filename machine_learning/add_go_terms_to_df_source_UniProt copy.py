import pandas as pd
import requests
from tqdm import tqdm
import os
import time
from urllib.parse import quote


# Fuzzy keyword matching
MITO_KEYWORDS = ["mitochondrion", "mitochondrial", "mito"]
ER_KEYWORDS = [
    "endoplasmic reticulum", "sarcoplasmic reticulum", "golgi", "lysosome",
    "granule", "membrane stack", "plasma membrane", "vacuole", "ER", "endoplasmic reticulum",
    "rough endoplasmic reticulum",
    "smooth endoplasmic reticulum",
    "sarcoplasmic reticulum",
    "er membrane",
    "er lumen",
    "er cisterna",
    "er tubule",
    "plasmodesmatal endoplasmic reticulum",
    "nuclear envelope",  # contiguous with ER
    "membrane stack",  # can refer to ER-like stacks
    "post-ER compartment",
    "secretory granule",
    "secretory vesicle",
    "pre-golgi compartment",
    "golgi intermediate compartment",
    "intermediate compartment",
    "perinuclear space",
    "endoplasmic reticulum-golgi intermediate compartment",  # ERGIC
    "perinuclear region",
    "coated vesicle",  # e.g., COPI/COPII
    "exit site",       # ER exit site
    "vesicle",
    "transport vesicle",
    "er exit site",
    "microsome"  # isolated ER vesicles
]

def classify_location(locations):
    locations = [loc.lower() for loc in locations]
    is_mito = any(any(k in loc for k in MITO_KEYWORDS) for loc in locations)
    is_er = any(any(k in loc for k in ER_KEYWORDS) for loc in locations)

    if is_mito and is_er:
        return "Multiple"
    elif is_mito:
        return "GO:0005739"
    elif is_er:
        return "GO:0005783"
    else:
        return "cyto_nuclear"

def fetch_uniprot_batch(accession_list):
    """
    Uses GET with accession:Pxxx OR accession:Qxxx ... for up to ~100 IDs per query.
    """
    base_url = "https://rest.uniprot.org/uniprotkb/search"
    headers = {"Accept": "application/json"}

    query_string = " OR ".join([f"accession:{acc}" for acc in accession_list])
    encoded_query = quote(query_string)

    url = f"{base_url}?query={encoded_query}&format=json&size={len(accession_list)}"

    response = requests.get(url, headers=headers)
    if response.status_code != 200:
        print(f"❌ Error {response.status_code}: {response.text}")
        return {}

    data = response.json()
    result = {}

    for entry in data.get("results", []):
        acc = entry["primaryAccession"]
        locs = []
        for comment in entry.get("comments", []):
            if comment.get("commentType") == "SUBCELLULAR LOCATION":
                for loc_entry in comment.get("subcellularLocations", []):
                    loc = loc_entry.get("location", {}).get("value")
                    if loc:
                        locs.append(loc)
        result[acc] = locs

    return result


# === Main Loop ===

general_working_dir = "pipeline/output/output_20250603_145910_ml_all_organisms"
organism_names = [
    "Homo_sapiens", "Mus_musculus", "Dario_rerio", "Daphnia_magna",
    "Caenorhabditis_elegans", "Drosophila_Melanogaster", "Arabidopsis_thaliana",
    "Physcomitrium_patens", "Chlamydomonas_reinhardtii",
    "Candida_glabrata", "Saccharomyces_cerevisiae", "Zygosaccharomyces_rouxii"
]

organism_names = ["Homo_sapiens"]  # For testing
batch_size = 400

for name in tqdm(organism_names, leave=False):
    input_file = os.path.join(general_working_dir, name, "feature_matrix_with_go_terms.csv")
    output_file = os.path.join(general_working_dir, name, "feature_matrix_with_go_terms_alt.csv")

    feature_matrix = pd.read_csv(input_file)
    feature_matrix1 = feature_matrix.copy()
    feature_matrix1['GO_Term'] = None

    ids = feature_matrix1['protein_id'].dropna().unique().tolist()
    go_mapping = {}
    # Batching with tqdm
    for i in tqdm(range(0, len(ids), batch_size), desc="Fetching UniProt data", leave=False):
        batch_ids = ids[i:i+batch_size]
        batch_result = fetch_uniprot_batch(batch_ids)
        go_mapping.update(batch_result)
        time.sleep(1)  # Respect UniProt rate limits

    # Apply classification to feature_matrix1
    feature_matrix1["GO_Term"] = feature_matrix1['protein_id'].apply(
        lambda pid: classify_location(go_mapping.get(pid, []))
    )

    # Save output
    feature_matrix1.to_csv(output_file, index=False)

    # Optional check
    matching_rows = feature_matrix1['GO_Term'] == feature_matrix.get('GO_Term', 'MISSING')
    percentage_match = (matching_rows.sum() / len(feature_matrix1)) * 100
    print(f"✅ {name}: {percentage_match:.2f}% GO_Term match (if column existed)")
