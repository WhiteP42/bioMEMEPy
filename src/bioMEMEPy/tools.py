from itertools import islice
import hashlib
import random


# Class PWM
class PWM:
    def __init__(self, seq: str, alphabet, m_length, top_val):
        self.alphabet = alphabet
        self.length = m_length
        self.matrix = {nucl: [float(0)] * self.length for nucl in alphabet}
        for i in range(self.length):
            for nucl in self.alphabet:
                if nucl == seq[i]:
                    self.matrix[nucl][i] = top_val
                else:
                    self.matrix[nucl][i] = (1 - top_val) / len(self.alphabet - 1)

    def normalize(self):
        for i in range(self.length):
            total = 0
            for nucl in self.alphabet:
                total += self.matrix[nucl][i]
            for nucl in self.alphabet:
                self.matrix[nucl][i] /= total

    def print(self):
        print(self.matrix)


# Class RPM
class RPM:
    def __init__(self, alphabet, m_length):
        self.alphabet = alphabet
        self.m_length = m_length
        self.matrix = dict()
        self.hash_map = dict()

    def add_seq(self, seq):
        hash_key = hashlib.sha256(seq.encode()).hexdigest()[:16]
        self.matrix[hash_key] = [float(0)] * (len(seq) - self.m_length + 1)
        self.hash_map[hash_key] = seq
        return hash_key

    def update_z(self, hash_key, z, offset):
        self.matrix[hash_key][offset] = z

    def normalize_seq(self, hash_key):
        total = 0
        for z in self.matrix[hash_key]:
            total += z
        for z in self.matrix[hash_key]:
            z /= total


def p0_gen(seqs, alphabet):
    prop = dict()
    total = 0
    for nucl in alphabet:
        prop[nucl] = 0
    for seq in seqs:
        for nucl in seq:
            prop[nucl] += 1
            total += 1
    for nucl in prop:
        prop[nucl] /= total
    return prop


# Load sequences from the FASTA file
def extract(filepath):
    seqs = []
    with open(filepath, 'r') as f:
        for line in islice(f, 1, None, 2):
            seqs.append(line.rstrip('\n'))
    return seqs


# Extract motif from sequence
def snip(seq, length, s_pos):
    if s_pos > len(seq) - s_pos + 1:
        raise ValueError('Snippet overflows the provided sequence.')

    f_pos = s_pos + length
    snippet = []
    for pos in range(s_pos, f_pos):
        snippet.append(seq[pos])
    return ''.join(snippet)


# Seeding functions
def gather(seqs, m_length, amount=0):
    ret_seqs = []

    if amount == 0:
        for seq in seqs:
            for pos in range(len(seq) - m_length):
                ret_seqs.append(snip(seq, m_length, pos))

    elif amount > 0:
        for n_gather in range(amount):
            seq = None
            while seq not in ret_seqs:
                seq = seqs[random.randint(0, len(seqs) - 1)]
            ret_seqs.append(snip(seq, m_length, random.randint(0, len(seq) - m_length)))

    else:
        raise ValueError('Amount must be 0 or higher.')
    return ret_seqs

# Get a hash
def get_hash(seq):
    return hashlib.sha256(seq.encode()).hexdigest()[:16]