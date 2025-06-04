import os
import subprocess
import MDAnalysis as mda
from MDAnalysis.analysis import secondary_structure
import matplotlib.pyplot as plt
import numpy as np

def run(cmd, cwd=None):
    print(f"Running: {cmd}")
    subprocess.run(cmd, shell=True, check=True, cwd=cwd)

def prepare_system():
    run("gmx pdb2gmx -f peptide_ace.pdb -o processed.gro -p topol.top -ff amber99sb-ildn -water tip3p")
    run("gmx editconf -f processed.gro -o boxed.gro -c -d 1.0 -bt cubic")
    run("gmx solvate -cp boxed.gro -cs spc216.gro -o solvated.gro -p topol.top")
    run("gmx grompp -f mdp/ions.mdp -c solvated.gro -p topol.top -o ions.tpr")
    run("echo 'SOL' | gmx genion -s ions.tpr -o ionized.gro -p topol.top -pname NA -nname CL -neutral")

def run_minimization():
    run("gmx grompp -f mdp/minim.mdp -c ionized.gro -p topol.top -o em.tpr")
    run("gmx mdrun -v -deffnm em")

def run_equilibration():
    run("gmx grompp -f mdp/nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr")
    run("gmx mdrun -deffnm nvt")
    run("gmx grompp -f mdp/npt.mdp -c nvt.gro -r nvt.gro -p topol.top -o npt.tpr")
    run("gmx mdrun -deffnm npt")

def run_production():
    run("gmx grompp -f mdp/md.mdp -c npt.gro -p topol.top -o md.tpr")
    run("gmx mdrun -deffnm md")

def analyze_helix(nterm_residues=12):
    u = mda.Universe("md.tpr", "md.xtc")
    dssp = secondary_structure.DSSP(u).run()
    ss = dssp.secondary_structure[:, :nterm_residues]
    helix_fraction = (ss == 'H').mean(axis=1)

    avg_fraction = helix_fraction.mean()
    return avg_fraction

def plot_helix_fraction():
    prepare_system()
    run_minimization()
    run_equilibration()
    run_production()
    analyze_helix(nterm_residues=12)

