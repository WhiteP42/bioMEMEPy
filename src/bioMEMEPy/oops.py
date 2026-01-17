import logging
import math
from . import tools
from .tools import BasePWM as PWM

logger = logging.getLogger(__name__)


class RPM(tools.BaseRPM):
    def softmax(self, hash_key):
        max_log = max(self.log_matrix[hash_key])
        for index, log in enumerate(self.log_matrix[hash_key]):
            self.resp_matrix[hash_key][index] = math.exp(log - max_log)
        total = sum(self.resp_matrix[hash_key])
        for index, exp in enumerate(self.resp_matrix[hash_key]):
            self.resp_matrix[hash_key][index] /= total


def e_step(pwm: PWM, rpm: RPM, seqs, p0):
    for seq in seqs:
        # Add a sequence and generate a hash for it on RPM.
        hash_key = rpm.add_seq(seq)
        # For a particular sequence, for every possible motif-length snip:
        for offset in range(len(seq) - pwm.length + 1):
            snippet = tools.snip(seq, pwm.length, offset)
            log = 0
            # For every position on the snip, calculate and accumulate log probability against background frequency.
            for index, nucl in enumerate(snippet):
                log += math.log(pwm.get_val(nucl, index)) - math.log(p0[nucl])
            rpm.update_log(hash_key, log, offset)
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
                for offset in range(len(seq) - pwm.length + 1):
                    snippet = tools.snip(seq, pwm.length, offset)
                    # If the nucleotide is the same as the target, add the responsibility value for that snippet.
                    if snippet[pos] == trgt_nucl:
                        count += rpm.resp_matrix[hash_key][offset]
            pwm.update(count, trgt_nucl, pos)
    # All columns should sum 1.
    pwm.normalize()

    # Rebuild background probabilities.
    new_p0 = {nucl: 0 for nucl in pwm.alphabet}
    z_val = {nucl: 0 for nucl in pwm.alphabet}
    total = tools.nucl_count(seqs, pwm.alphabet)
    for seq in seqs:
        hash_key = tools.get_hash(seq)
        for offset in range(len(seq) - pwm.length + 1):
            snippet = tools.snip(seq, pwm.length, offset)
            for i, nucl in enumerate(snippet):
                z_val[nucl] += rpm.resp_matrix[hash_key][offset]
    for nucl in pwm.alphabet:
        new_p0[nucl] = total[nucl] - z_val[nucl]
    p0.clear()
    p0.update(new_p0)


def oops(seqs, alphabet, m_length, top_val, extract_val, threshold, max_iter):
    logger.debug('Called model OOPS.')
    # Get required k-mers to seed.
    logger.debug('Gathering seed candidates.')
    if tools.snip_count(seqs, m_length) >= 10000:
        seed_seqs = tools.gather(seqs, m_length, extract_val)
    else:
        seed_seqs = tools.gather(seqs, m_length)
    logger.debug('Done!')
    logger.debug(f'Total seeds: {len(seed_seqs)}.')

    # Seeding process:
    top_candidate = None
    logger.debug(f'Beginning seeding process.')
    for index, snip in enumerate(seed_seqs):
        logger.debug(f'Testing {snip} ({index + 1}/{len(seed_seqs)}).')
        p0 = tools.p0_gen(seqs, alphabet)
        current_pwm = PWM(snip, alphabet, m_length, top_val)
        current_rpm = RPM(m_length)
        # Run 1 EM interation with the seed.
        e_step(current_pwm, current_rpm, seqs, p0)
        m_step(current_pwm, current_rpm, seqs, p0)

        # Compute total log-likelihood and compare with the previous best (or generate best).
        log_like = tools.log_like(current_rpm)
        if top_candidate is None or log_like > top_candidate[1]:
            logger.debug('Replacing top candidate.')
            top_candidate = [snip, log_like]
        logger.debug(f'Log-like computed: {log_like}. Top candidate is {top_candidate[1]}.')

    # Run EM to convergence with the selected seed and return PWM.
    logger.debug('Beginning EM...')
    seed = top_candidate[0]
    p0 = tools.p0_gen(seqs, alphabet)
    pwm = PWM(seed, alphabet, m_length, top_val)
    rpm = RPM(m_length)
    converged = False
    prev_loglike = None
    iterations = 0
    while not converged:
        iterations += 1
        logger.debug(f'Iteration {iterations}.')
        e_step(pwm, rpm, seqs, p0)
        m_step(pwm, rpm, seqs, p0)
        log_like = tools.log_like(rpm)
        if prev_loglike is None or (log_like - prev_loglike) > threshold:
            prev_loglike = log_like
        elif iterations >= max_iter or (log_like - prev_loglike) <= threshold:
            converged = True

    return pwm.matrix