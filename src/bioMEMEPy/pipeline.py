import logging
from .oops import oops
from .tools import extract

def meme(fasta, alphabet, motif_l, model):
    # Errors
    if not isinstance(motif_l, int):
        raise TypeError('Length of the motif must be int.')

    #Read FASTA
    with open(fasta, 'r') as f:
        seqs = extract(fasta)

    #Call mode
    _models = {'oops': oops} # Future models to be implemented here.
    try:
        result = _models[model](seqs, alphabet, motif_l)
    except KeyError:
        raise ValueError(f'Model {model} unsupported.')