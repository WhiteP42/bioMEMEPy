import logging
import math
from . import tools
from. tools import BasePWM as PWM


class RPM(tools.BaseRPM):
    def softmax(self, hash_key):
        max_log = max(self.matrix[hash_key])
        for log in self.matrix[hash_key]:
            log = math.exp(log - max_log)

    def normalize_seq(self, hash_key):
        total = 0
        for z in self.matrix[hash_key]:
            total += z
        for z in self.matrix[hash_key]:
            z /= total


def e_step(pwm: PWM, rpm: RPM, p0, seqs):
    for seq in seqs:
        hash_key = rpm.add_seq(seq)
        for offset in range(len(seq) - pwm.length + 1):
            snippet = tools.snip(seq, pwm.length, offset)
            for j in range(len(snippet)):
                nucl = snippet[j]
                log_nucl = math.log(pwm.matrix[nucl][j]) - math.log(p0[nucl])
            rpm.update_val(hash_key, log_nucl, offset)
        rpm.softmax(hash_key)
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