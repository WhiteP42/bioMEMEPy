from itertools import islice

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
        self.matrix = {nucl: [0] * length for nucl in alphabet}