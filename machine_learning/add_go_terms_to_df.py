import os
import pandas as pd
from goatools.obo_parser import GODag
from goatools.base import download_go_basic_obo


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

def load_obo():
    '''
    load the latest obo file to exchange the GO terms with their aspects
    '''
    obo_path = "pipeline/go.obo"
    if not os.path.exists(obo_path):
        obo_path = download_go_basic_obo()  
    # Load the GO DAG
    go_dag = GODag(obo_path)
    return go_dag

def run(name):
    # Load the GO ontology file (e.g., gene_ontology.obo)
    go_obo = load_obo()

    go_dag = GODag("go-basic.obo")  


    go_ids = [
        #Go terms for: ER, Golgi, Ribosome, Mitochondria, Nucleus, Lysosome, Cell membrane, Cytoplasm
        "GO:0005783",
        "GO:0005794",
        "GO:0005840",
        "GO:0005739",
        "GO:0005634",
        "GO:0005764",
        "GO:0005886",
        "GO:0005737"
    ]
    # Set the working directory
    working_dir = os.path.dirname("pipeline/output/output_20250519_142700_machine_learning_human/" + name + "/")
    feature_matrix_path = working_dir + "/feature_matrix.csv"
    # Read the feature matrix from the CSV file
    feature_matrix = pd.read_csv(feature_matrix_path, index_col=0)
    print(feature_matrix.head())
    # Add GO terms to the feature matrix

    go_child_dir = working_dir + "/go_childs/"
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

    # Add filtered GO terms rowwise to the feature matrix
    for protein_id, row in feature_matrix.iterrows():
        filtered_go_annotations = {}
        filtered_terms = []
        if protein_id in go_annotations:
            for term in go_annotations[protein_id]:
                if term in child_dict:
                    filtered_terms.append(child_dict[term])
            feature_matrix.at[protein_id, "GO_Term"] = ",".join(filtered_terms)


    '''
    for protein_id, row in feature_matrix.iterrows():
        filtered_terms = []
        if protein_id in go_annotations:
            for term in go_annotations[protein_id]:
                # test if term has any children terms
                if term in go_dag:
                    children = go_dag[term].children
                    if not children:
                        filtered_terms.append(term)
            feature_matrix.at[protein_id, "GO_Term"] = ",".join(filtered_terms)
    '''
    # Filter rows to keep only those containing the GO term GO:0005739
    #feature_matrix = feature_matrix[feature_matrix["GO_Term"].str.contains("GO:0005739", na=False)]

    # Remove duplicate GO terms in each row
    feature_matrix["GO_Term"] = feature_matrix["GO_Term"].apply(
        lambda x: ",".join(sorted(set(x.split(",")))) if pd.notna(x) else x
    )
    '''
    # Remove the GO term GO:0005739 from the GO_Term column
    feature_matrix["GO_Term"] = feature_matrix["GO_Term"].apply(
        lambda x: ",".join(term for term in x.split(",") if term != "GO:0005739") if pd.notna(x) else x
    )'''
    '''
    # if a protein has multiple go terms, copy the row and assign the go term to the new row
    # and remove the go term from the original row
    feature_matrix["GO_Term"] = feature_matrix["GO_Term"].apply(
        lambda x: x.split(",") if pd.notna(x) else x
    )'''

    feature_matrix = feature_matrix.explode("GO_Term")
    # Empty the cells that have more than one GO term
    feature_matrix["GO_Term"] = feature_matrix["GO_Term"].apply(
        lambda x: "" if isinstance(x, str) and "," in x else x)
    # Remove rows with an empty GO term
    '''for index, row in feature_matrix.iterrows():
        if isinstance(row["GO_Term"], str):
            terms = row["GO_Term"].split(",")
            if len(terms) == 1 and terms[0] == "":
                feature_matrix.drop(index, inplace=True)'''

    # moves the go term column to the beginning of the DataFrame
    go_term_column = feature_matrix.pop("GO_Term")
    feature_matrix.insert(0, "GO_Term", go_term_column)
    # Save the updated feature matrix to a CSV file
    output_path = working_dir + "/feature_matrix_with_go_terms.csv"
    feature_matrix.to_csv(output_path)


