from goatools.obo_parser import GODag
import os


"""
This script retrieves all child terms of cellular anatomical structures from the GO DAG.
It includes terms that are part of the cellular component namespace and checks for 'part_of' relationships
to ensure comprehensive coverage of related terms.
It also handles the 'occurs_in' relationships and intersection_of entries that may indicate part_of relationships.
The results are saved in a specified directory, with each term's children listed in a separate file
"""


# Load the GO DAG (Gene Ontology terms)
go_dag = GODag("pipeline/go.obo")  # Adjust the path to the correct OBO file
def get_part_of_children(go_id, go_dag, visited=None):
    if visited is None:
        visited = set()
    children = []
    for term in go_dag.values():
        # Check is_a relationships
        if hasattr(term, 'is_a') and go_id in term.is_a:
            if term.id not in visited:
                children.append(term.id)
                visited.add(term.id)
                children.extend(get_part_of_children(term.id, go_dag, visited))
        # Check relationships dictionary for part_of and occurs_in
        if hasattr(term, 'relationships'):
            for rel_type, rel_list in term.relationships.items():
                if rel_type in ['part_of', 'occurs_in']:
                    if go_id in rel_list:
                        if term.id not in visited:
                            children.append(term.id)
                            visited.add(term.id)
                            children.extend(get_part_of_children(term.id, go_dag, visited))
        # Check intersection_of for part_of relationships
        if hasattr(term, 'intersection_of'):
            for entry in term.intersection_of:
                # entry can be a tuple or a string
                if isinstance(entry, tuple):
                    if entry[0] == 'part_of' and entry[1] == go_id:
                        if term.id not in visited:
                            children.append(term.id)
                            visited.add(term.id)
                            children.extend(get_part_of_children(term.id, go_dag, visited))
                elif isinstance(entry, str):
                    # Sometimes the entry is a string like "part_of GO:0005739"
                    if entry.strip().startswith('part_of') and go_id in entry:
                        if term.id not in visited:
                            children.append(term.id)
                            visited.add(term.id)
                            children.extend(get_part_of_children(term.id, go_dag, visited))
    return children


def get_all_children(term):
    children = []
    def _get_all_children(t):
        for child in t.children:
            children.append(child.id)
            _get_all_children(child)
    _get_all_children(term)
    return children
def get_all_descendants(go_id, go_dag):
        children = []
        for term in go_dag.values():
            if term.namespace == "cellular_component":
                if go_id in term.get_all_parents():
                    children.append(term.id)
        return children
go_ids = [
    #Go terms for: ER, Golgi, Ribosome, Mitochondria, Nucleus, Lysosome, Cell membrane, Cytoplasm
    "GO:0005783",
    "GO:0005794",
    "GO:0005840",
    "GO:0005739",
    "GO:0005634",
    "GO:0005764",
    "GO:0005886",
    "GO:0005737",
    "GO:0005829",
]

for go_id in go_ids:
    
    # Get the term from the GO DAG
    go_term = go_dag[go_id]
    print(f"Cellular anatomical structure term: {go_term}")

    # Print out the children of the cellular anatomical structure term
    for child in go_term.children:
        print(f"Child term: {child.id} - {child.name}")

    # Recursively get all child terms and part_of children
    all_children = get_all_children(go_term)
    part_of_children = get_part_of_children(go_id, go_dag)
    more_children = get_all_descendants(go_id, go_dag)

    # Remove duplicates by converting to a set
    child_terms = set(all_children + part_of_children + more_children)

    # deltete duplicates
    child_terms = set(child_terms)


    # Print out the number of child terms: {len(child_terms)}")

    # Print out the child terms
    for go_id_child in child_terms:
        term = go_dag[go_id_child]
        print(f"{go_id_child}: {term.name}")

    # Create a directory for the child terms
    os.makedirs("pipeline/output/output_20250519_142700_machine_learning_human/Homo_sapiens/go_childs", exist_ok=True)
    # Save the sorted child_terms
    with open(f"pipeline/output/output_20250519_142700_machine_learning_human/Homo_sapiens/go_childs/{go_id}.txt", "w") as file:
        if not child_terms:
            file.write(f"No children found for {go_id} - {go_term.name}\n")
        for go_id_child in sorted(child_terms):
            term = go_dag[go_id_child]
            file.write(f"{go_id_child}_{term.name}\n")
        file.write(f"{go_id}_{go_term.name}")
        if not child_terms:
            file.write(f"No children found for {go_id} - {go_term.name}\n")