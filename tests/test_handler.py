from bioMEMEPy.handler import PWM
from bioMEMEPy.handler import snip
from bioMEMEPy import monomers as mnm

def test_pwm():
    pwm = PWM(mnm.dna, 4)
    expected_empty_pwm = {'A': [0, 0, 0, 0],
                          'T': [0, 0, 0, 0],
                          'C': [0, 0, 0, 0],
                          'G': [0, 0, 0, 0]}
    assert pwm.matrix == expected_empty_pwm
    pwm.init()
    for i in range(pwm.length):
        total = 0
        for nucl in mnm.dna:
            total += pwm.matrix[nucl][i]
        assert total == 1 # Handle the float error!!

def snip_test():
    seq = 'ATACGTTAT'
    snippet = snip(seq, 4, 4)
    assert snippet == 'GTTA'