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
    reference_name = max(sequences, key=lambda x: len(sequences[x]))
    reference_seq = sequences[reference_name]
    aligned_sequences = {}

    for name, seq in sequences.items():
        if name == reference_name:
            aligned_sequences[name] = reference_seq
        else:
            # Use scoring scheme: match=2, mismatch=-1, open_gap=-0.5, extend_gap=-0.1
            alignment = pairwise2.align.globalms(reference_seq, seq, 2, -1, -0.5, -0.1, one_alignment_only=True)[0]
            aligned_sequences[name] = alignment[1]  # aligned target sequence
    return aligned_sequences


# Plot the aligned profiles as a heatmap
def plot_property_profiles(aligned_sequences):
    label_map = {"A": 0, "B": 1, "N": 2, "H": 3, "X": 4, "-": np.nan}
    cmap = plt.get_cmap("tab10", 5)

    max_length = min(100, max(len(seq) for seq in aligned_sequences.values()))  # Limit to first 100 amino acids
    numeric_profiles = {}

    for name, seq in aligned_sequences.items():
        numeric_seq = [label_map.get(c, np.nan) for c in seq[:100]]  # Slice to first 100 amino acids
        padded_seq = numeric_seq + [np.nan] * (max_length - len(numeric_seq))
        numeric_profiles[name] = padded_seq

    _, ax = plt.subplots(figsize=(14, len(numeric_profiles)))  # Removed unused variable 'fig'
    for i, (name, profile) in enumerate(numeric_profiles.items()):
        profile_array = np.array(profile, dtype=float)
        ax.imshow([profile_array], aspect="auto", cmap=cmap, extent=[0, max_length, i, i + 1])

    ax.set_yticks(np.arange(len(numeric_profiles)) + 0.5)
    ax.set_yticklabels(numeric_profiles.keys())
    ax.set_xlabel("Aligned Position")
    ax.set_title("Aligned Biochemical Property Profiles (First 100 Amino Acids)")

    sm = plt.cm.ScalarMappable(cmap=cmap)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, ticks=range(5))
    cbar.ax.set_yticklabels(["Acidic", "Basic", "Neutral", "Hydrophobic", "Other"])
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
    print(f"Aligned sequences saved to {output_path}")    
plt.show()  # Show the heatmap plot