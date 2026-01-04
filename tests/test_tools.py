import pytest
from bioMEMEPy.tools import snip
from bioMEMEPy import mnm

def test_snip():
    seq = 'ATACGTTAT'
    assert snip(seq, 5, 1) == 'TACGT'
    assert snip(seq, 4, 5) == 'TTAT'
    with pytest.raises(ValueError):
        snip(seq, 4, 6)