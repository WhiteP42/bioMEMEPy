from itertools import islice
import random

# Load sequences from the FASTA file
def extract(filepath):
    seqs = []
    with open(filepath, 'r') as f:
        for line in islice(f, 1, None, 2):
            seqs.append(line.rstrip('\n'))
    return seqs

# Class PWM
class PWM:
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
                self.matrix[nucl][i] = random.random()
        self.normalize()

    def update(self, seq):
        for i in range(self.length):
            self.matrix[seq[i]][i] += 1