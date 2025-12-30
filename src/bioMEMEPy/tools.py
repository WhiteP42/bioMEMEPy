from itertools import islice
from random import random

# Load sequences from the FASTA file
def extract(filepath):
    seqs = []
    with open(filepath, 'r') as f:
        for line in islice(f, 1, None, 2):
            seqs.append(line.rstrip('\n'))
    return seqs

# Extract motif from sequence
def snip(seq, length, s_pos):
    if s_pos > len(seq) - s_pos:
        raise ValueError('Snippet overflows the provided sequence.')

    f_pos = s_pos + length
    snippet = []
    for pos in range(s_pos, f_pos):
        snippet.append(seq[pos])
    return ''.join(snippet)

# Class PSSM
class PSSM:
    def __init__(self, alphabet, length):
        self.alphabet = alphabet
        self.length = length
        self.matrix = {nucl: [float(0)] * length for nucl in alphabet}

    def normalize(self):
        for i in range(self.length):
            total = 0
            for nucl in self.alphabet:
                total += self.matrix[nucl][i]
            for nucl in self.alphabet:
                self.matrix[nucl][i] /= total

    def init(self):
        for nucl in self.matrix:
            for i in range(self.length):
                self.matrix[nucl][i] = random()
        self.normalize()

    def print(self):
        print(self.matrix)