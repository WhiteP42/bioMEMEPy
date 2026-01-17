import math
from itertools import islice
import hashlib
import random
import logging

logger = logging.getLogger(__name__)


# Total count of each nucleotide:
def nucl_count(seqs, alphabet) -> dict:
    total = {nucl: 0 for nucl in alphabet}
    for nucl in alphabet:
        for seq in seqs:
            total[nucl] += seq.count(nucl)
    return  total


def p0_gen(seqs, alphabet):
    prop = nucl_count(seqs, alphabet)
    total = sum(prop.values())
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
    if s_pos > (len(seq) - length + 1):
        raise ValueError('Snippet overflows the provided sequence.')

    f_pos = s_pos + length
    snippet = []
    for pos in range(s_pos, f_pos):
        snippet.append(seq[pos])
    return ''.join(snippet)


# How many possible snippets in a sequence
def snip_count(seqs, m_length):
    count = 0
    for seq in seqs:
        count += (len(seq) - m_length + 1)
    return count


# Seeding functions
def gather(seqs, m_length, amount=0):
    ret_snips = []

    if amount == 0:
        for seq in seqs:
            logger.debug(f'Running sequence {seq}.')
            for pos in range(len(seq) - m_length + 1):
                logger.debug(f'Position: {pos}/{len(seq) - m_length}')
                ret_snips.append(snip(seq, m_length, pos))

    elif amount > 0:
        while len(ret_snips) < amount:
            seq = seqs[random.randint(0, len(seqs) - 1)]
            snippet = snip(seq, m_length, random.randint(0, len(seq) - m_length + 1))
            if snippet not in ret_snips:
                ret_snips.append(snippet)

    else:
        raise ValueError('Amount must be 0 or higher.')
    return ret_snips

# Get a hash
def get_hash(seq):
    return hashlib.sha256(seq.encode()).hexdigest()[:16]

# Compute log-likelihood:
def log_like(rpm):
    max_log = 0
    total_z = 0
    for log_vctr in rpm.log_matrix:
        for log_val in rpm.log_matrix[log_vctr]:
            if log_val > max_log:
                max_log = log_val
    for log_vctr in rpm.log_matrix:
        for log_val in rpm.log_matrix[log_vctr]:
            total_z += math.exp(log_val - max_log)
    return max_log + math.log(total_z)

#Obtain consensus from PWM:
def consensus(pwm: dict, m_length, alphabet):
    consens = []
    for index in range(m_length):
        best = None
        for nucl in alphabet:
            if best is None or pwm[nucl][index] > best[1]:
                best = [nucl, pwm[nucl][index]]
        consens.append(best[0])
    return ''.join(consens)


# Class PWM
class BasePWM:
    def __init__(self, seq: str, alphabet, m_length, top_val):
        self.alphabet = alphabet
        self.length = m_length
        self.matrix = {nucl: [float(0)] * self.length for nucl in alphabet}
        for i in range(self.length):
            for nucl in self.alphabet:
                if nucl == seq[i]:
                    self.matrix[nucl][i] = top_val
                else:
                    self.matrix[nucl][i] = (1 - top_val) / (len(self.alphabet) - 1)

    def update(self, value, nucl, pos):
        self.matrix[nucl][pos] = value

    def normalize(self):
        for i in range(self.length):
            total = 0
            for nucl in self.alphabet:
                total += self.matrix[nucl][i]
            for nucl in self.alphabet:
                self.matrix[nucl][i] /= total

    def get_val(self, nucl, pos):
        return self.matrix[nucl][pos]

    def print(self):
        print(self.matrix)


# Class RPM
class BaseRPM:
    def __init__(self, m_length):
        self.m_length = m_length
        self.resp_matrix = dict()
        self.log_matrix = dict()
        self.hash_map = dict()

    def add_seq(self, seq):
        hash_key = get_hash(seq)
        self.resp_matrix[hash_key] = [float(0)] * (len(seq) - self.m_length + 1)
        self.log_matrix[hash_key] = [float(0)] * (len(seq) - self.m_length + 1)
        self.hash_map[hash_key] = seq
        return hash_key

    def update_log(self, hash_key, val, offset):
        self.log_matrix[hash_key][offset] = val