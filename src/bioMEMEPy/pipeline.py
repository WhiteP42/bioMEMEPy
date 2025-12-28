import logging
from Bio import SeqIO
from .oops import oops

def meme(fasta, alphabet, motif_l, model):
    # Errors
    if not isinstance(motif_l, int):
        raise TypeError('Length of the motif must be int.')

    #Read FASTA
    with open(fasta, 'r') as f:
        seqs = list(SeqIO.parse(f, 'fasta'))

    #Call mode
    _models = {'oops': oops} # Future models to be implemented here.
    try:
        result = _models[model](seqs, alphabet, motif_l)
    except KeyError:
        raise ValueError(f'Model {model} unsupported.')