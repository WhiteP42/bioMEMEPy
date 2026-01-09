import logging
import math
from . import tools


def e_step(pwm: tools.PWM, rpm: tools.RPM, p0, seqs):
    for seq in seqs:
        hash_key = rpm.add_seq(seq)
        for offset in range(len(seq) - pwm.length + 1):
            z = 0
            snippet = tools.snip(seq, pwm.length, offset)
            for j in range(len(snippet)):
                nucl = snippet[j]
                log_nucl = math.log(pwm.matrix[nucl][j]) - math.log(p0[nucl])
                z += log_nucl
            rpm.update_z(hash_key, z, offset)
        rpm.normalize_seq(hash_key)


def m_step():
    raise NotImplementedError


def oops(seqs, alphabet, m_length, top_val=0.5, extract_val=2000):
    # Get all k-mers to seed.
    if len(seqs) * (len(seqs[0] - m_length + 1)) >= 10000:
        seed_seqs = tools.gather(seqs, m_length, extract_val)
    else:
        seed_seqs = tools.gather(seqs, m_length)

    # Begin seeding.
    for seq in seed_seqs:
        p0 = tools.p0_gen(seed_seqs, alphabet)
        current_pwm = tools.PWM(seq, alphabet, m_length, top_val)
        current_rpm = tools.RPM(alphabet, m_length)
        e_step(current_pwm, current_rpm, p0, seqs)