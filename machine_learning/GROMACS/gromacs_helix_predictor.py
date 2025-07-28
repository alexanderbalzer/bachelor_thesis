import os
import subprocess
import tempfile
from pathlib import Path

def run_helix_propensity(sequence: str) -> float:
    sequence = sequence.strip().upper()
    acetylated = sequence.startswith("ACE-")
    seq = sequence.replace("ACE-", "")

    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        pdb_path = tmpdir / "peptide.pdb"

        # Build initial structure using Biopython + PeptideBuilder
        from Bio.PDB import PDBIO
        from PeptideBuilder import Geometry, PeptideBuilder
        import Bio.PDB

        geo = Geometry.geometry(seq[0])
        structure = PeptideBuilder.initialize_res(geo)
        for aa in seq[1:]:
            PeptideBuilder.add_residue(structure, aa)
        io = PDBIO()
        io.set_structure(structure)
        with open(pdb_path, "w") as handle:
            io.save(handle)
        '''
        if acetylated:
            with open(pdb_path, "r") as f:
                lines = f.readlines()

            # Manually add ACE residue using standard coordinates (dummy, editable)
            ace_residue = [
                "HETATM    1  CH3 ACE     1      -0.500   0.000   0.000  1.00  0.00           C\n",
                "HETATM    2  C   ACE     1       0.000   0.000   0.000  1.00  0.00           C\n",
                "HETATM    3  O   ACE     1       1.200   0.000   0.000  1.00  0.00           O\n"
            ]

            lines = ace_residue + lines
            with open(pdb_path, "w") as f:
                f.writelines(lines)
        '''
        # Generate topology
        result = subprocess.run([
            "gmx", "pdb2gmx", "-f", str(pdb_path), "-o", "processed.gro", "-p", "topol.top", "-i", "posre.itp",
            "-ff", "amber99sb-ildn", "-water", "tip3p"
        ], input=b"1\n1\n", cwd=tmpdir, capture_output=True)
        if result.returncode != 0:
            print(result.stderr.decode())
            raise RuntimeError("GROMACS command failed")

        # Solvate, add ions, minimize, run MD, and analyze as before
        result = subprocess.run(["gmx", "editconf", "-f", "processed.gro", "-o", "newbox.gro", "-c", "-d", "1.0", "-bt", "cubic"], cwd=tmpdir, input=b"\n", capture_output=True)
        if result.returncode != 0:
            print(result.stderr.decode())
            raise RuntimeError("GROMACS command failed")

        result = subprocess.run(["gmx", "solvate", "-cp", "newbox.gro", "-cs", "spc216.gro", "-o", "solv.gro", "-p", "topol.top"], cwd=tmpdir, capture_output=True)
        if result.returncode != 0:
            print(result.stderr.decode())
            raise RuntimeError("GROMACS command failed")

        with open(tmpdir / "ions.mdp", "w") as f:
            f.write("integrator = steep\n")
        result = subprocess.run(["gmx", "grompp", "-f", "ions.mdp", "-c", "solv.gro", "-p", "topol.top", "-o", "ions.tpr"], cwd=tmpdir, capture_output=True)
        if result.returncode != 0:
            print(result.stderr.decode())
            raise RuntimeError("GROMACS command failed")

        result = subprocess.run(["gmx", "genion", "-s", "ions.tpr", "-o", "solv_ions.gro", "-p", "topol.top", "-pname", "NA", "-nname", "CL", "-neutral"], cwd=tmpdir, input=b"13\n", capture_output=True)
        if result.returncode != 0:
            print(result.stderr.decode())
            raise RuntimeError("GROMACS command failed")

        with open(tmpdir / "em.mdp", "w") as f:
            f.write("integrator = steep\nemtol = 1000.0\nnsteps = 50000\n")
        result = subprocess.run(["gmx", "grompp", "-f", "em.mdp", "-c", "solv_ions.gro", "-p", "topol.top", "-o", "em.tpr"], cwd=tmpdir, capture_output=True)
        if result.returncode != 0:
            print(result.stderr.decode())
            raise RuntimeError("GROMACS command failed")

        result = subprocess.run(["gmx", "mdrun", "-deffnm", "em"], cwd=tmpdir, capture_output=True)
        if result.returncode != 0:
            print(result.stderr.decode())
            raise RuntimeError("GROMACS command failed")

        with open(tmpdir / "md.mdp", "w") as f:
            f.write("""integrator = md
nsteps = 50000
dt = 0.002
nstxout = 500
nstvout = 500
nstenergy = 500
nstlog = 500
continuation = no
constraint_algorithm = lincs
constraints = all-bonds
cutoff-scheme = Verlet
ns_type = grid
rlist = 1.0
rcoulomb = 1.0
rvdw = 1.0
tcoupl = V-rescale
tc-grps = System
tau_t = 0.1
ref_t = 300
pcoupl = no
""")
        result = subprocess.run(["gmx", "grompp", "-f", "md.mdp", "-c", "em.gro", "-p", "topol.top", "-o", "md.tpr"], cwd=tmpdir, capture_output=True)
        if result.returncode != 0:
            print(result.stderr.decode())
            raise RuntimeError("GROMACS command failed")

        result = subprocess.run(["gmx", "mdrun", "-deffnm", "md"], cwd=tmpdir, capture_output=True)
        if result.returncode != 0:
            print(result.stderr.decode())
            raise RuntimeError("GROMACS command failed")
        
        from Bio.PDB import PDBParser, DSSP

        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("peptide", pdb_path)
        model = structure[0]

        dssp = DSSP(model, pdb_path)
        helix_residues = [res_id for res_id in dssp if dssp[res_id][2] == 'H']
        helix_fraction = len(helix_residues) / len(dssp)
        
        return helix_fraction * 100


        xpm_path = tmpdir / "ss.xpm"
        if not xpm_path.exists():
            raise RuntimeError("DSSP analysis failed")

        with open(xpm_path) as f:
            data = f.read()
        helix_count = data.count('H')
        total_count = sum(data.count(c) for c in ['H', 'E', 'C', 'T', 'S', 'B', 'G', 'I'])
        if total_count == 0:
            raise ValueError("No secondary structure found")
        return round(100 * helix_count / total_count, 2)


if __name__ == "__main__":
    # Example usage
    sequence = "ACE-MKLLLILLLLILLLL"
    try:
        helix_percentage = run_helix_propensity(sequence)
        print(f"Helix percentage: {helix_percentage}%")
    except Exception as e:
        print(f"Error: {e}")