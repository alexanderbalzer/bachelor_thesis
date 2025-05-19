from goatools.obo_parser import GODag

# Load the GO DAG (Gene Ontology terms)
go_dag = GODag("go-basic.obo")  # Adjust the path to the correct OBO file

# Define the GO term for 'cellular anatomical structure'
cellular_anatomical_structure_id = "GO:0008104"

# Get the term from the GO DAG
cellular_anatomical_structure_term = go_dag[cellular_anatomical_structure_id]
print(f"Cellular anatomical structure term: {cellular_anatomical_structure_term}")
# children
# Print out the children of the cellular anatomical structure term
for child in cellular_anatomical_structure_term.children:
     print(f"Child term: {child.id} - {child.name}")
# Retrieve all child terms of this GO term
#child_terms = cellular_anatomical_structure_term.get_children()


# Initialize a list to store direct child terms
direct_children = []

for child in cellular_anatomical_structure_term.children:
    # Check if the child term is a direct child (not a grandchild)
    if child.id.startswith("GO:"):
        direct_children.append(child.id)

child_terms = direct_children

# Print out the child terms
for go_id in child_terms:
    term = go_dag[go_id]
    print(f"{go_id}: {term.name}")
# save the sorted child_terms
with open("pipeline/output/output_20250519_142700_machine_learning_human/Homo_sapiens/targeting_anchoring_targeting_transport.txt", "w") as file:
    for go_id in sorted(child_terms):
        term = go_dag[go_id]
        file.write(f"{go_id}_{term.name}\n")