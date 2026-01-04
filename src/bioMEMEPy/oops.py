import logging
from . import tools


def e_step(pwm: tools.PWM, rpm: tools.RPM, seqs):
    for seq in seqs:
        rpm.add_seq(seq)


def m_step():
    raise NotImplementedError


def oops(seqs, alphabet, m_length, top_val=0.5, extract_val=2000):
    # Get all k-mers to seed.
    if len(seqs) * (len(seqs[0] - m_length)) >= 10000:
        seed_seqs = tools.gather(seqs, m_length, extract_val)
    else:
        seed_seqs = tools.gather(seqs, m_length)

    # Begin seeding.
    for seq in seed_seqs:
        current_pwm = tools.PWM(seq, alphabet, m_length, top_val)
        current_rpm = tools.RPM(alphabet, m_length)
        e_step(current_pwm, current_rpm, seqs)