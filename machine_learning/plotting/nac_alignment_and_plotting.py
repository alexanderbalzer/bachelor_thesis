import matplotlib.pyplot as plt
import numpy as np
from Bio import SeqIO
# Removed unused import: SequenceMatcher
import os
from Bio import SeqIO, pairwise2
# Removed unused import: format_alignment
import subprocess


def simplify_sequence_to_properties(seq):
    acidic = set("ED")
    basic = set("KRH")
    neutral = set("GPSTA")
    hydrophobic = set("LIVFMYWC")

    simplified = []
    for aa in seq:
        if aa in acidic:
            simplified.append("A")
        elif aa in basic:
            simplified.append("B")
        elif aa in neutral:
            simplified.append("N")
        elif aa in hydrophobic:
            simplified.append("H")
        else:
            simplified.append("X")
    return "".join(simplified)

def load_fasta(fasta_path):
    sequences = {}
    for record in SeqIO.parse(fasta_path, "fasta"):
        sequences[record.id] = str(record.seq)
    return sequences

def load_and_simplify_fasta(fasta_path):
    sequences = {}
    for record in SeqIO.parse(fasta_path, "fasta"):
        simplified = simplify_sequence_to_properties(str(record.seq))
        sequences[record.id] = simplified
    return sequences

# Step 3: Align all sequences to the longest one
def align_to_reference(sequences):
    reference_name = 'sp|E9PAV3|NACAM_HUMAN'
    reference_seq = sequences[reference_name]
    aligned_sequences = {}

    for name, seq in sequences.items():
        if name == reference_name:
            aligned_sequences[name] = reference_seq
        else:
            # Use scoring scheme: match=2, mismatch=-1, open_gap=-0.5, extend_gap=-0.2
            alignment = pairwise2.align.globalms(reference_seq, seq, 9, 1, -8, -7, one_alignment_only=True)[0]
            aligned_seq = alignment[1]  # aligned target sequence

            # Replace mismatches with '.'
            aligned_seq = "".join(c if c != "-" else "." for c in aligned_seq)
            aligned_sequences[name] = aligned_seq  #

    return aligned_sequences


# Plot the aligned profiles as a heatmap
def plot_property_profiles(aligned_sequences):
    import pandas as pd

    label_map = {"A": 0, "B": 1, "N": 2, "H": 3, "X": 4, "-": np.nan, ".": np.nan, " ": np.nan, "\t": np.nan}
    cmap = plt.get_cmap("tab10", 5)

    max_length = max(len(seq) for seq in aligned_sequences.values())  # Use full length

    # Build DataFrame: rows=sequence names, columns=aligned positions
    data = {}
    for name, seq in aligned_sequences.items():
        numeric_seq = [label_map.get(c, np.nan) for c in seq[:max_length]]
        # Pad to max_length
        if len(numeric_seq) < max_length:
            numeric_seq += [np.nan] * (max_length - len(numeric_seq))
        data[name] = numeric_seq
    df = pd.DataFrame.from_dict(data, orient='index', columns=np.arange(1, max_length+1))

    fig, ax = plt.subplots(figsize=(14, len(df)))
    im = ax.imshow(df.values, aspect="auto", cmap=cmap, interpolation='none', vmin=0, vmax=4)

    ax.set_yticks(np.arange(len(df)))
    ax.set_yticklabels(df.index)
    # Show every 10th position, but cover the full length
    step = max(1, max_length // 20)
    ax.set_xticks(np.arange(0, max_length, step))
    ax.set_xticklabels(np.arange(1, max_length+1, step))
    ax.set_xlabel("Aligned Position")
    ax.set_title("Aligned Biochemical Property Profiles (Full Sequence)")

    sm = plt.cm.ScalarMappable(cmap=cmap)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax)
    cbar.set_ticklabels(["Acidic (A)", "Basic (B)", "Neutral (N)", "Hydrophobic (H)", "Unknown (X)"])
    cbar.set_label("Residue Type")

    plt.tight_layout()
    
def run_iupred(sequence, tmp_fasta="tmp_iupred.fasta"):
    iupred = '/home/abalzer/iupred2a/iupred2a.py'
    with open(tmp_fasta, "w") as f:
        f.write(">query\n" + sequence + "\n")
    result = subprocess.run(["python3", iupred, tmp_fasta, "short"],
                            stdout=subprocess.PIPE, text=True)
    disorder_scores = []
    for line in result.stdout.splitlines():
        if line.strip().startswith("#") or len(line.strip().split()) < 3:
            continue
        try:
            _, _, score = line.strip().split()[:3]  # Removed unused variables 'idx' and 'aa'
            disorder_scores.append(float(score))
        except ValueError:
            continue
    return disorder_scores
def find_morf_candidates(sequence, disorder_scores, min_length=8, threshold=0.5):  # Removed unused parameter 'max_length'
    candidates = []
    start = None
    for i, score in enumerate(disorder_scores):
        if score > threshold:
            if start is None:
                start = i
        else:
            if start is not None and (i - start) >= min_length:
                candidates.append((start, i))  # (start, end)
            start = None
    # check end
    if start is not None and (len(disorder_scores) - start) >= min_length:
        candidates.append((start, len(disorder_scores)))
    return [simplify_sequence_to_properties(sequence[start:end]) for (start, end) in candidates]

def compare_morfs(morfs_dict):
    from difflib import SequenceMatcher
    consensus = []
    keys = list(morfs_dict.keys())
    for i, seq1 in enumerate(morfs_dict[keys[0]]):
        for seq2 in morfs_dict[keys[1]]:
            if SequenceMatcher(None, seq1, seq2).ratio() > 0.7:
                consensus.append((seq1, seq2))
    return consensus


# Example usage
if __name__ == "__main__":
    working_dir = "pipeline/output/output_20250617_183139_latest_ML/NACA isoforms and homologes"  
    fasta_path = os.path.join(working_dir, "naca_isoforms_and_homologues.fasta")
    seqs = load_fasta(fasta_path)
    morfs = {}
    # Snip the last 192 amino acids off each sequence
    for name in seqs:
        seqs[name] = seqs[name][:-192]
    for name, full_seq in seqs.items():
        disorder_scores = run_iupred(full_seq)
        morfs[name] = find_morf_candidates(full_seq, disorder_scores)
    conserved = compare_morfs(morfs)
    print("Conserved MoRFs:")
    for name, seq in morfs.items():
        print(f"{name}: {', '.join(seq)}")
        # Save conserved MoRFs to a file
        conserved_morfs_path = os.path.join(working_dir, "conserved_morfs.txt")
        with open(conserved_morfs_path, "w") as conserved_file:
            for name, seq in morfs.items():
                conserved_file.write(f"{name}: {', '.join(seq)}\n")
        print(f"Conserved MoRFs saved to {conserved_morfs_path}")
    print("Conserved MoRF pairs:")
    for a, b in conserved:
        print(f"Conserved MoRF: {a} ~ {b}")
    sequences = load_and_simplify_fasta(fasta_path)
    aligned = align_to_reference(sequences)
    plot_property_profiles(aligned)
    # Save the aligned sequences to a file
    output_path = os.path.join(working_dir, "aligned_sequences.fasta")
    with open(output_path, "w") as output_file:
        for name, seq in aligned.items():
            output_file.write(f">{name}\n{seq}\n")
    plt.savefig(os.path.join(working_dir, "aligned_property_profiles.pdf"), dpi=500, bbox_inches="tight")
    print(f"Aligned sequences saved to {output_path}")