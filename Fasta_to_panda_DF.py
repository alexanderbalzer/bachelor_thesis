import fastaparser as fp
import matplotlib as plt
import numpy as np
import pandas as pd


def fasta_to_dataframe(fasta_file):
    data = []
    with open(fasta_file, 'r') as file:
        for sequence in fp.Reader(file):
            data.append([sequence.id, sequence.sequence])
        return pd.DataFrame(data, columns=['id', 'sequence'])

print(fasta_to_dataframe('P53_HUMAN.fasta'))