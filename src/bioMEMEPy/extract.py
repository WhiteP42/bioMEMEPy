from itertools import islice

# Load sequences from the FASTA file

def extract(file):
    seq = []
    with open(file, 'r') as f:
        for line in islice(f, 0, None, 2):
            seq.append(line.rstrip('\n'))
    return seq