import pytest
from bioMEMEPy import tools
from bioMEMEPy import mnm

def test_snip():
    seq = 'ATACGTTAT'
    assert tools.snip(seq, 5, 1) == 'TACGT'
    assert tools.snip(seq, 4, 5) == 'TTAT'
    with pytest.raises(ValueError):
        tools.snip(seq, 4, 6)

def test_gather_0():
    seqs = ['AAAA', 'TTTT', 'CCCC', 'GGGG']
    assert tools.gather(seqs, 3) == [
        'AAA', 'AAA',
        'TTT', 'TTT',
        'CCC', 'CCC',
        'GGG', 'GGG'
    ]

def test_gather_3():
    seqs = ['AAAA', 'TTTT', 'CCCC', 'GGGG']
    selected = tools.gather(seqs, 3, 3)
    assert len(selected) == 3