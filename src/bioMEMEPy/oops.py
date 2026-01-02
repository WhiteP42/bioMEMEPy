import logging
from . import tools

def oops(seqs, alphabet, m_length, top_val=0.5, extract_val=2000):
    # Get all k-mers to seed.
    if len(seqs) * (len(seqs[0] - m_length)) >= 10000:
        seed_seqs = tools.gather(seqs, m_length, extract_val)
    else:
        seed_seqs = tools.gather(seqs, m_length)