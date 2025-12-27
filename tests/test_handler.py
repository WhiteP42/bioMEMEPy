from bioMEMEPy.tools import PSSM
from bioMEMEPy.tools import snip
from bioMEMEPy import monomers as mnm

def test_pssm():
    pssm = PSSM(mnm.dna, 4)
    expected_empty_pssm = {'A': [0, 0, 0, 0],
                          'T': [0, 0, 0, 0],
                          'C': [0, 0, 0, 0],
                          'G': [0, 0, 0, 0]}
    assert pssm.matrix == expected_empty_pssm
    pssm.init()
    for i in range(pssm.length):
        total = 0
        for nucl in mnm.dna:
            total += pssm.matrix[nucl][i]
        assert total == 1 # Handle the float error!!

def test_snip():
    seq = 'ATACGTTAT'
    snippet = snip(seq, 4, 4)
    assert snippet == 'GTTA' # Wrong format!!