import logging
import time
from .oops import oops
from .tools import extract

# Logger call:
logger = logging.getLogger(__name__)


def meme(fasta, alphabet, model, m_length, motif_num=1, top_val=0.5, extract_val=2000, threshold=1e-5, max_iter=200):
    # Start time
    start_time = time.perf_counter()

    # Errors
    if not isinstance(m_length, int):
        raise TypeError('Length of the motif must be int.')

    #Read FASTA
    with open(fasta, 'r') as f:
        seqs = extract(f)

    #Call mode
    _models = {'oops': oops} # ZOOPS and ANR models to be added here.
    try:
        result = _models[model](seqs, alphabet, m_length, top_val, extract_val, threshold, max_iter)
    except KeyError:
        raise ValueError(f'Model {model} unsupported.')

    # Runtime
    total_time = time.perf_counter() - start_time
    logger.info(f'Runtime: %.3f s.', total_time)

    #Return PWM
    return result