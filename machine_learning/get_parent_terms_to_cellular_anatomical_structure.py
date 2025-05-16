from goatools.obo_parser import GODag

# Load the GO DAG (Gene Ontology terms)
go_dag = GODag("go-basic.obo")  # Adjust the path to the correct OBO file

# Define the GO term for 'cellular anatomical structure'
cellular_anatomical_structure_id = "GO:0110165"

# Get the term from the GO DAG
cellular_anatomical_structure_term = go_dag[cellular_anatomical_structure_id]

# Retrieve all child terms of this GO term
child_terms = cellular_anatomical_structure_term.get_all_children()

# Print out the child terms
for go_id in sorted(child_terms):
    term = go_dag[go_id]
    print(f"{go_id}: {term.name}")
# save the sorted child_terms
with open("pipeline/output/output_20250515_105213/Homo_sapiens/cellular_child_terms.txt", "w") as file:
    for go_id in sorted(child_terms):
        term = go_dag[go_id]
        file.write(f"{go_id}\n")