#!/usr/bin/env python

"""
Calculates a set of properties from a protein sequence:
    - hydrophobicity (according to a particular scale)
    - mean hydrophobic dipole moment assuming it is an alpha-helix.
    - total charge (at pH 7.4)
    - amino acid composition
    - discimination factor according to Rob Keller (IJMS, 2011)

Essentially the same as HeliQuest (reproduces the same values).

Author:
  Joao Rodrigues
  j.p.g.l.m.rodrigues@gmail.com
"""

from __future__ import print_function

import argparse
import csv
import math
import os
import time
import numpy as np
import pandas as pd

#
# Definitions
#
scales = {'Fauchere-Pliska': {'A':  0.31, 'R': -1.01, 'N': -0.60,
                              'D': -0.77, 'C':  1.54, 'Q': -0.22,
                              'E': -0.64, 'G':  0.00, 'H':  0.13,
                              'I':  1.80, 'L':  1.70, 'K': -0.99,
                              'M':  1.23, 'F':  1.79, 'P':  0.72,
                              'S': -0.04, 'T':  0.26, 'W':  2.25,
                              'Y':  0.96, 'V':  1.22},

          'Eisenberg': {
    'A':  0.62,  # Alanine
    'R': -2.53,  # Arginine
    'N': -0.78,  # Asparagine
    'D': -0.90,  # Aspartic acid
    'C':  0.29,  # Cysteine
    'Q': -0.85,  # Glutamine
    'E': -0.74,  # Glutamic acid
    'G':  0.48,  # Glycine
    'H': -0.40,  # Histidine
    'I':  1.38,  # Isoleucine
    'L':  1.06,  # Leucine
    'K': -1.50,  # Lysine
    'M':  0.64,  # Methionine
    'F':  1.19,  # Phenylalanine
    'P':  0.12,  # Proline
    'S': -0.18,  # Serine
    'T': -0.05,  # Threonine
    'W':  0.81,  # Tryptophan
    'Y':  0.26,  # Tyrosine
    'V':  1.08   # Valine
},
          }
_supported_scales = list(scales.keys())

aa_charge = {'E': -1, 'D': -1, 'K': 1, 'R': 1}

#
# Functions
#
def assign_hydrophobicity(sequence, scale='Eisenberg'):  # noqa: E302
    """Assigns a hydrophobicity value to each amino acid in the sequence"""

    hscale = scales.get(scale, None)
    if not hscale:
        raise KeyError('{} is not a supported scale. '.format(scale))

    hvalues = []
    for aa in sequence:
        sc_hydrophobicity = hscale.get(aa, None)
        if sc_hydrophobicity is None:
            raise KeyError('Amino acid not defined in scale: {}'.format(aa))
        hvalues.append(sc_hydrophobicity)

    return hvalues


def calculate_moment(array, angle=100):
    """Calculates the hydrophobic dipole moment from an array of hydrophobicity
    values. Formula defined by Eisenberg, 1982 (Nature). Returns the average
    moment (normalized by sequence length)

    uH = sqrt(sum(Hi cos(i*d))**2 + sum(Hi sin(i*d))**2),
    where i is the amino acid index and d (delta) is an angular value in
    degrees (100 for alpha-helix, 180 for beta-sheet).
    """

    sum_cos, sum_sin = 0.0, 0.0
    for i, hv in enumerate(array):
        rad_inc = ((i*angle)*math.pi)/180.0
        sum_cos += hv * math.cos(rad_inc)
        sum_sin += hv * math.sin(rad_inc)
    return math.sqrt(sum_cos**2 + sum_sin**2)


def calculate_charge(sequence, charge_dict=aa_charge):
    """Calculates the charge of the peptide sequence at pH 7.4
    """
    sc_charges = [charge_dict.get(aa, 0) for aa in sequence]
    return sum(sc_charges)

def assign_charge(sequence, acetylated, seq_range):
    """Assigns charge values to each amino acid in the sequence"""
    charges = []
    for i, aa in enumerate(sequence):
        if i == 0 and seq_range == 0 and not acetylated:
            charges.append(1)  # if the first amino acid is not acetylated and NatA/D, it has a positive charge
        elif i == 0 and seq_range == 0 and acetylated:
            charges.append(0)
        elif i == 0 and seq_range == 0 and len(sequence) > 1 and sequence[1] in {'D', 'E', 'N', 'K', 'L', 'I', 'F', 'W'} and not acetylated:
            charges.append(1)
        elif i == 0 and seq_range == 0 and len(sequence) > 1 and sequence[1] in {'D', 'E', 'N', 'K', 'L', 'I', 'F', 'W'} and acetylated:
            charges.append(0)
        elif aa == 'H':
            charges.append(0.1)  # Approximate partial charge at pH 7.4
        else:
            charges.append(np.abs(aa_charge.get(aa, 0)))
    return charges


def calculate_vector_moment(values, angle=100):
    """Generic vector moment calculator (used for hydrophobic or electrostatic)"""
    sum_cos, sum_sin = 0.0, 0.0
    for i, val in enumerate(values):
        rad_inc = ((i * angle) * math.pi) / 180.0
        sum_cos += val * math.cos(rad_inc)
        sum_sin += val * math.sin(rad_inc)
    return sum_cos, sum_sin


def calculate_alignment(hvec, qvec):
    """Compute cosine similarity between hydrophobic and electrostatic vectors"""
    hmag = math.sqrt(hvec[0]**2 + hvec[1]**2)
    qmag = math.sqrt(qvec[0]**2 + qvec[1]**2)
    if hmag == 0 or qmag == 0:
        return 0.0  # Avoid division by zero
    dot = hvec[0]*qvec[0] + hvec[1]*qvec[1]
    return dot / (hmag * qmag)



def calculate_discrimination(mean_uH, total_charge):
    """Returns a discrimination factor according to Rob Keller (IJMS, 2011)
    A sequence with d>0.68 can be considered a potential lipid-binding region.
    """
    d = 0.944*mean_uH + 0.33*total_charge
    return d


def calculate_composition(sequence):
    """Returns a dictionary with percentages per classes"""

    # Residue character table
    polar_aa = set(('S', 'T', 'N', 'H', 'Q', 'G'))
    speci_aa = set(('P', 'C'))
    apolar_aa = set(('A', 'L', 'V', 'I', 'M'))
    charged_aa = set(('E', 'D', 'K', 'R'))
    aromatic_aa = set(('W', 'Y', 'F'))

    n_p, n_s, n_a, n_ar, n_c = 0, 0, 0, 0, 0
    for aa in sequence:
        if aa in polar_aa:
            n_p += 1
        elif aa in speci_aa:
            n_s += 1
        elif aa in apolar_aa:
            n_a += 1
        elif aa in charged_aa:
            n_c += 1
        elif aa in aromatic_aa:
            n_ar += 1

    return {'polar': n_p, 'special': n_s,
            'apolar': n_a, 'charged': n_c, 'aromatic': n_ar}

helix_propensity = {
    'A': 1.45, 'C': 0.77, 'D': 0.98, 'E': 1.53,
    'F': 1.12, 'G': 0.53, 'H': 1.24, 'I': 1.00,
    'K': 1.07, 'L': 1.34, 'M': 1.20, 'N': 0.73,
    'P': 0.59, 'Q': 1.17, 'R': 0.79, 'S': 0.79,
    'T': 0.82, 'V': 1.14, 'W': 1.14, 'Y': 0.61
}

def helix_scoring(seq, helix_propensity=helix_propensity):
    helix_score = sum(helix_propensity.get(res.upper(), 0) for res in seq) / len(seq)
    return helix_score

def analyze_sequence(name=None, sequence=None, window=9, verbose=False, w_h = 0.944 / (0.944 + 0.33), w_q = 0.33 / (0.944 + 0.33)):
    """Runs all the above on a sequence. Pretty prints the results"""

    if not sequence:
        raise Exception('Either I need glasses or there is no sequence.')

    if not name:
        name = 'Unnamed'
    if sequence.startswith('X'):
        acetylated = True
        sequence = sequence[1:]  # remove the acetylation prefix
    else:
        acetylated = False
    w = window
    if w < 0:
        w = len(sequence)  # automatically set

    outdata = []  # for csv writing

    # Processing...
    seq_len = len(sequence)
    #print('[+] Analysing sequence {} ({} aa.)'.format(name, seq_len))
    #print('[+] Using a window of {} aa.'.format(w))
    for seq_range in range(1, 10):
        if seq_range + w > seq_len:
            if verbose:
                print(f'Skipping window {seq_range+1} due to insufficient length: {seq_len} < {seq_range+w}')
            break
        seq_w = sequence[seq_range:seq_range+w]
        if seq_range and len(seq_w) < w:
            break
        
        helix_score = helix_scoring(seq_w, helix_propensity)
        if helix_score < 0.7:
            if verbose:
                print(f'Skipping window {seq_range+1} due to low helix score: {helix_score:.2f}')
            _t = [name, sequence, seq_range+1, w, seq_w, 0, 0, 0, 0, 0,
                  0, 0, 0, 0, 0, 0, helix_score]
            outdata.append(_t)
            continue

        # Numerical values
        z = calculate_charge(seq_w)
        seq_h = assign_hydrophobicity(seq_w)
        av_h = sum(seq_h)/len(seq_h)
        
        h_cos, h_sin = calculate_vector_moment(seq_h)
        av_uH = math.sqrt(h_cos**2 + h_sin**2)/len(seq_h)
        
        seq_q = assign_charge(seq_w, acetylated, seq_range)
        seq_q = - np.abs(np.array(seq_q))
        q_cos, q_sin = calculate_vector_moment(seq_q)
        av_uQ = math.sqrt(q_cos**2 + q_sin**2)/len(seq_q)

        # Calculate the keller weighted linear combination of the two vectors
        combined_magnitude = w_h * av_uH + w_q * av_uQ  
        alignment = calculate_alignment((h_cos, h_sin), (q_cos, q_sin))
        #alignment = np.abs(alignment)  # Ensure non-negative value

        d = calculate_discrimination(av_uH, z)


        # AA composition
        aa_comp = calculate_composition(seq_w)
        n_tot_pol = aa_comp['polar'] + aa_comp['charged']
        n_tot_apol = aa_comp['apolar'] + aa_comp['aromatic'] + aa_comp['special']  # noqa: E501
        n_charged = aa_comp['charged']  # noqa: E501
        n_aromatic = aa_comp['aromatic']  # noqa: E501

        _t = [name, sequence, seq_range+1, w, seq_w, z, av_h, av_uH, av_uQ, alignment, d,
      n_tot_pol, n_tot_apol, n_charged, n_aromatic, combined_magnitude, helix_score]

        outdata.append(_t)

        if verbose:
            print('  Window {}: {}-{}-{}'.format(seq_range+1, seq_range,
                                                 seq_w, seq_range+w))
            print('    z={:<3d} <H>={:4.3f} <uH>={:4.3f} D={:4.3f}'.format(z, av_h,  # noqa: E501
                                                                           av_uH, d))  # noqa: E501
            print('    Amino acid composition')
            print('      Polar    : {:3d} / {:3.2f}%'.format(n_tot_pol, n_tot_pol*100/w))  # noqa: E501
            print('      Non-Polar: {:3d} / {:3.2f}%'.format(n_tot_apol, n_tot_apol*100/w))  # noqa: E501
            print('      Charged  : {:3d} / {:3.2f}%'.format(n_charged, n_charged*100/w))  # noqa: E501
            print('      Aromatic : {:3d} / {:3.2f}%'.format(n_aromatic, n_aromatic*100/w))  # noqa: E501
            print()

    return outdata

def analyze_sequence_with_set_parameters(name=None, sequence=None, seq_range=0, w=18, verbose=False, w_h = 0.944 / (0.944 + 0.33), w_q = 0.33 / (0.944 + 0.33)):

    if sequence.startswith('X'):
        acetylated = True
        sequence = sequence[1:]  # remove the acetylation prefix
    else:
        acetylated = False
    seq_w = sequence[seq_range:seq_range+w]

    helix_score = helix_scoring(seq_w, helix_propensity)
    # Numerical values
    z = calculate_charge(seq_w)
    seq_h = assign_hydrophobicity(seq_w)
    av_h = sum(seq_h)/len(seq_h)
    
    h_cos, h_sin = calculate_vector_moment(seq_h)
    av_uH = math.sqrt(h_cos**2 + h_sin**2)/len(seq_h)
    
    seq_q = assign_charge(seq_w, acetylated, seq_range)
    seq_q = - np.abs(np.array(seq_q))
    q_cos, q_sin = calculate_vector_moment(seq_q)
    av_uQ = math.sqrt(q_cos**2 + q_sin**2)/len(seq_q)

    # Calculate the keller weighted linear combination of the two vectors
    combined_magnitude = w_h * av_uH + w_q * av_uQ  
    alignment = calculate_alignment((h_cos, h_sin), (q_cos, q_sin))

    d = calculate_discrimination(av_uH, z)


    # AA composition
    aa_comp = calculate_composition(seq_w)
    n_tot_pol = aa_comp['polar'] + aa_comp['charged']
    n_tot_apol = aa_comp['apolar'] + aa_comp['aromatic'] + aa_comp['special']  # noqa: E501
    n_charged = aa_comp['charged']  # noqa: E501
    n_aromatic = aa_comp['aromatic']  # noqa: E501

    _t = [name, sequence, seq_range+1, w, seq_w, z, av_h, av_uH, av_uQ, alignment, d,
    n_tot_pol, n_tot_apol, n_charged, n_aromatic, combined_magnitude, helix_score]

    return _t

def read_fasta_file(afile):
    """Parses a file with FASTA formatted sequences"""

    if not os.path.isfile(afile):
        raise IOError('File not found/readable: {}'.format(afile))

    sequences = []
    seq_name, cur_seq = None, None
    with open(afile) as handle:
        for line in handle:
            line = line.strip()
            if line.startswith('>'):
                if cur_seq:
                    sequences.append((seq_name, ''.join(cur_seq)))
                seq_name = line[1:]
                cur_seq = []
            elif line:
                cur_seq.append(line)
    sequences.append((seq_name, ''.join(cur_seq)))  # last seq

    return sequences

def run(sequence, verbose):
    data = []
    w_h = 0.944 #/ (0.944 + 0.33) # weight for hydrophobic vector
    w_q = 0.33 #/ (0.944 + 0.33) # weight for charge vector
    for i in range(10, 120):
        temp_data = analyze_sequence(sequence=sequence, verbose=verbose, window=i, w_h=w_h, w_q=w_q)
        data.extend(temp_data)
    # _t = [name, sequence, seq_range+1, w, seq_w, z, av_h, av_uH, av_uQ, alignment, d,
    #  n_tot_pol, n_tot_apol, n_charged, n_aromatic, combined_magnitude]
    av_uH = [row[7] for row in data]  # This is the average hydrophobic moment
    if av_uH is None or len(av_uH) == 0:
        max_av_uH = 0
        start_best_window = 0
        length_best_window = 0
        electrostatic_help = 0
        discrimination_factor = 0
        helix_score = 0
        charge = 0
    else:
        max_av_uH = max(av_uH)
        max_index = av_uH.index(max_av_uH)
        start_best_window = data[max_index][2]  # This is the start index of the best window
        length_best_window = data[max_index][3]  # This is the length of the best window
        av_uQ = data[max_index][8] # This is the average electrostatic moment
        alignment = data[max_index][9]  # This is the alignment value
        electrostatic_help = alignment * av_uQ  # This is the electrostatic help
        discrimination_factor = data[max_index][10]
        helix_score = data[max_index][16]  # This is the helix score
        charge = data[max_index][5]  # This is the charge of the best window
    df = pd.DataFrame(data, columns=['Name', 'Sequence', 'Start', 'Window Size', 'Sub-Sequence',
                                       'Charge', 'Mean Hydrophobicity',
                                       'Mean Hydrophobic Moment', 'Mean Electrostatic Moment',
                                       'Alignment', 'Discrimination Factor',
                                       'No. Polar AA', 'No. Apolar AA',
                                       'No. Charged AA', 'No. Aromatic AA',
                                       'Combined Magnitude', 'Helix Score'])
    '''with open('hydrophobicity_results.csv', 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(df.columns)
        writer.writerows(df.values)'''
    return max_av_uH, start_best_window, length_best_window, electrostatic_help, discrimination_factor, helix_score, charge

def run_alternative(sequence, name='Unnamed', seq_range=0, w=18, verbose=False):
    data = analyze_sequence_with_set_parameters(name=name, sequence=sequence, seq_range=seq_range, w=w, verbose=verbose)
    av_uQ = data[8] # This is the average electrostatic moment
    alignment = data[9]  # This is the alignment value
    electrostatic_help = alignment * av_uQ  # This is the electrostatic help
    return electrostatic_help
'''_t = [name, sequence, seq_range+1, w, seq_w, z, av_h, av_uH, av_uQ, alignment, d,
    n_tot_pol, n_tot_apol, n_charged, n_aromatic, combined_magnitude, helix_score]'''
if __name__ == '__main__':
    print(run("MLR", verbose=False))  # Example sequence


'''
    ap = argparse.ArgumentParser(description=__doc__)
    i_opts = ap.add_mutually_exclusive_group(required=True)
    i_opts.add_argument('-s', '--sequence',
                        help='Sequence of amino acids to analyze')
    i_opts.add_argument('-f', '--seqfile',
                        help='File with sequences in FASTA format')
    ap.add_argument('-o', '--outfile',
                    help='File to write results in CSV format')
    ap.add_argument('-v', '--verbose', action='store_true',
                    help='Write information to screen as well')
    ap.add_argument('--scale', choices=_supported_scales,
                    help='Hydrophobicity scale to use')
    ap.add_argument('-w', '--window', default=18, type=int,
                    help=(
                        'AA window to use during analysis. Set to -1 to '
                        'automatically match the full-length of the sequence'))
    cmd = ap.parse_args()

    # File or Sequence?
    if cmd.sequence:
        all_data = analyze_sequence(sequence=cmd.sequence, window=cmd.window,
                                    verbose=cmd.verbose)
    elif cmd.seqfile:
        seq_list = read_fasta_file(cmd.seqfile)
        all_data = []
        for name, seq in seq_list:
            data = analyze_sequence(name=name, sequence=seq, window=cmd.window,
                                    verbose=cmd.verbose)
            all_data += data'''


'''
    if not cmd.outfile:
        if cmd.seqfile:
            root, _ = os.path.splitext(cmd.seqfile)
        else:
            root = 'seq'
        outfn = root + '_hydrophobicity.txt'
    else:
        outfn = cmd.outfile

    print('[+] Writing results to {}'.format(outfn))
    if os.path.isfile(outfn):
        root, _ = os.path.splitext(outfn)
        outfn = root + '_' + time.strftime('%H%M%S%d%m%Y') + '.txt'
        print('  File already exists')
        print('  Writing to {}'.format(outfn))

    with open(outfn, 'w') as csvfile:
        _header = ['Name', 'Sequence', 'Window', 'Window Size', 'Sub-Sequence',
                   'Charge', 'Mean Hydrophobicity', 'Mean Hydrophobic Moment',
                   'Discrimination Factor', 'No. Polar AA', 'No. Apolar AA',
                   'No. Charged AA', 'No. Aromatic AA']

        writer = csv.writer(csvfile, dialect='excel')
        writer.writerow(_header)
        writer.writerows(all_data)
'''