from bioMEMEPy.tools import snip
from bioMEMEPy import mnm

def test_snip():
    seq = 'ATACGTTAT'
    snippet = snip(seq, 4, 4)
    assert snippet == 'GTTA'