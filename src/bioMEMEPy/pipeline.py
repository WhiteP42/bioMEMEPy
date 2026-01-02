import logging
import time
from .oops import oops
from .tools import extract

def meme(fasta, alphabet, motif_l, model):
    # Start time and logger config
    logger = logging.getLogger(__name__)
    start_time = time.perf_counter()

    # Errors
    if not isinstance(motif_l, int):
        raise TypeError('Length of the motif must be int.')

    #Read FASTA
    with open(fasta, 'r') as f:
        seqs = extract(f)

    #Call mode
    _models = {'oops': oops} # ZOOPS and ANR models to be added here.
    try:
        result = _models[model](seqs, alphabet, motif_l)
    except KeyError:
        raise ValueError(f'Model {model} unsupported.')

    # Runtime
    total_time = time.perf_counter() - start_time
    logger.info(f'Runtime: %.3f s.', total_time)