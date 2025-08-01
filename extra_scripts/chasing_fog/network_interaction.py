import requests
import networkx as nx
from pyvis.network import Network
from tqdm import tqdm
import time

# Seed proteins from NAC and TOM complexes
seed_proteins = {
    "NAC": ["NACA", "BTF3"],
    "TOM": ["TOMM20", "TOMM22", "TOMM40", "TOMM70A"]
}

# UniProt interaction endpoint (STRING is another option)
def get_interactors(uniprot_id):
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"
    try:
        response = requests.get(url)
        if response.status_code == 200:
            data = response.json()
            interactors = []
            if "comments" in data:
                for c in data["comments"]:
                    if c["commentType"] == "INTERACTION":
                        for interact in c.get("interactions", []):
                            interactor = interact["interactantTwo"]["uniProtKBAccession"]
                            interactors.append(interactor)
            return interactors
        else:
            return []
    except:
        return []

# Recursive interaction fetcher
def fetch_interactions(seeds, depth=1, delay=1.0):
    graph = nx.Graph()
    visited = set()

    def recurse(current_ids, d):
        if d == 0:
            return
        next_ids = set()
        for pid in tqdm(current_ids, desc=f"Depth {depth-d+1}"):
            if pid in visited:
                continue
            visited.add(pid)
            interactors = get_interactors(pid)
            for inter in interactors:
                graph.add_edge(pid, inter)
                if inter not in visited:
                    next_ids.add(inter)
            time.sleep(delay)
        recurse(next_ids, d-1)

    recurse(seeds, depth)
    return graph

# Collect all unique UniProt IDs (you can query UniProt to get these from gene names)
seed_ids = [*seed_proteins["NAC"], *seed_proteins["TOM"]]

# Fetch interaction network
interaction_graph = fetch_interactions(seed_ids, depth=2)

# Visualization
net = Network(notebook=True, height="750px", width="100%", bgcolor="#222222", font_color="white")
net.from_nx(interaction_graph)

# Highlight original complexes
for n in seed_proteins["NAC"]:
    if n in net.nodes:
        net.get_node(n)["color"] = "blue"
for n in seed_proteins["TOM"]:
    if n in net.nodes:
        net.get_node(n)["color"] = "red"

# Save and show
net.show("nac_tom_network.html")
