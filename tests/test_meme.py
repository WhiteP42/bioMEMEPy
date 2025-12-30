from bioMEMEPy import meme
from bioMEMEPy import mnm

def test_meme_oops():
    result = meme('q1.fasta', mnm.dna, 10, 'oops')