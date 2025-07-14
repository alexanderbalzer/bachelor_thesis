from Bio.PDB import MMCIFParser, PDBIO, PDBParser
import os


os.chdir("/home/abalzer/Downloads/fold_nacad_htt")
cif_file = "pdb3ior.ent"
pdb_file = "pdb3ior.pdb"

cif_parser = MMCIFParser(QUIET=True)
parser = PDBParser()
structure = parser.get_structure("3IOR", cif_file)

io = PDBIO()
io.set_structure(structure)
io.save(pdb_file)
print(f"Saved PDB to {pdb_file}")
