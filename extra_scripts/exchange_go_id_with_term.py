import os
import re
from goatools.obo_parser import GODag

# Path to your parent directory containing GO ID folders
parent_dir = "pipeline/output/output_20250515_105213/Homo_sapiens/go_term_rf_results"
# Path to your OBO file
obo_path = "pipeline/go.obo"

# Load GO DAG
go_dag = GODag(obo_path)

# Helper to get aspect letter from GO ID
def get_go_aspect(go_id, go_dag):
    """
    Get the aspect of a GO term.
    """
    if go_id in go_dag:
        go_term = go_dag[go_id]
        return go_term.name
    else:
        return "Unknown"

def sanitize_filename(name):
    # Replace any character that is not alphanumeric or underscore with underscore
    return re.sub(r'[^A-Za-z0-9_]', '_', name)

# Iterate through folders and rename
for folder in os.listdir(parent_dir):
    folder_path = os.path.join(parent_dir, folder)
    if os.path.isdir(folder_path) and folder.startswith("C_GO:"):
        go_id = folder[2:]  # Remove the "C_" prefix
        aspect = get_go_aspect(go_id, go_dag)
        aspect = aspect.replace(" ", "_")
        aspect = sanitize_filename(aspect)
        new_folder_name = f"{aspect}_{go_id}"
        new_folder_path = os.path.join(parent_dir, new_folder_name)
        os.rename(folder_path, new_folder_path)
        print(f"Renamed {folder} -> {new_folder_name}")
    else:
        print(f"Skipping {folder} (not a valid GO term folder)")
        continue

# After renaming, delete empty folders
for folder in os.listdir(parent_dir):
    folder_path = os.path.join(parent_dir, folder)
    if os.path.isdir(folder_path) and not os.listdir(folder_path):
        os.rmdir(folder_path)
        print(f"Deleted empty folder: {folder}")
