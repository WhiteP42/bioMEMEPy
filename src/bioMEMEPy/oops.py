import logging
from .pipeline import logger
import math
from . import tools
from .tools import BasePWM as PWM


class RPM(tools.BaseRPM):
    def softmax(self, hash_key):
        max_log = max(self.matrix[hash_key])
        for index, log in enumerate(self.matrix[hash_key]):
            self.matrix[hash_key][index] = math.exp(log - max_log)
        total = sum(self.matrix[hash_key])
        for index, log in enumerate(self.matrix[hash_key]):
            self.matrix[hash_key][index] /= total


def e_step(pwm: PWM, rpm: RPM, seqs, p0):
    for seq in seqs:
        # Add a sequence and generate a hash for it on RPM.
        hash_key = rpm.add_seq(seq)
        # For a particular sequence, for every possible motif-length snip:
        for offset in range(tools.snip_count(seq, pwm.length)):
            snippet = tools.snip(seq, pwm.length, offset)
            log = 0
            # For every position on the snip, calculate and accumulate log probability against background frequency.
            for index, nucl in enumerate(snippet):
                log += math.log(pwm.get_val(nucl, index)) - math.log(p0[nucl])
            rpm.update(hash_key, log, offset)
        rpm.softmax(hash_key)


def m_step(pwm: PWM, rpm: RPM, seqs, p0):
    # Update PWM with weighted counts.
    # For each position in the motif:
    for pos in range(pwm.length):
        # And for every nucleotide:
        for trgt_nucl in pwm.alphabet:
            count = 0
            # Go through every sequence and every snippet:
            for seq in seqs:
                hash_key = tools.get_hash(seq)
                for offset in range(tools.snip_count(seq, pwm.length)):
                    snippet = tools.snip(seq, pwm.length, offset)
                    # If the nucleotide is the same as the target, add the responsibility value for that snippet.
                    if snippet[pos] == trgt_nucl:
                        count += rpm.matrix[hash_key][offset]
            pwm.update(count, trgt_nucl, pos)
    # All columns should sum 1.
    pwm.normalize()

    # Rebuild background probabilities.
    new_p0 = {nucl: 0 for nucl in pwm.alphabet}
    z_val = {nucl: 0 for nucl in pwm.alphabet}
    total = tools.nucl_count(seqs, pwm.alphabet)
    for seq in seqs:
        hash_key = tools.get_hash(seq)
        for offset in range(tools.snip_count(seq, pwm.length)):
            snippet = tools.snip(seq, pwm.length, offset)
            for i, nucl in enumerate(snippet):
                z_val[nucl] += rpm.matrix[hash_key][offset]
    for nucl in pwm.alphabet:
        new_p0[nucl] = total[nucl] - z_val[nucl]
    p0.clear()
    p0.update(new_p0)


def oops(seqs, alphabet, m_length, top_val=0.5, extract_val=2000):
    # Get required k-mers to seed.
    if tools.snip_count(seqs, m_length) >= 10000:
        seed_seqs = tools.gather(seqs, m_length, extract_val)
    else:
        seed_seqs = tools.gather(seqs, m_length)

    # Seeding process:
    for seq in seed_seqs:
        # Generate background frequencies.
        p0 = tools.p0_gen(seed_seqs, alphabet)
        # Generate PWM and responsibility matrix (RPM).
        current_pwm = PWM(seq, alphabet, m_length, top_val)
        current_rpm = RPM(m_length)
        # Do 1 EM loop.
        e_step(current_pwm, current_rpm, seqs, p0)
        m_step(current_pwm, current_rpm, seqs, p0)
        # TODO: Add checking, save if better seed than previous or no seed is saved.
        # TODO: Run EM to convergence with the selected seed.