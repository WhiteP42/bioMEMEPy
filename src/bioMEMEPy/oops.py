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
        self.log_like += max_log + math.log(sum(self.resp_matrix[hash_key]))
        shift_sum = sum(self.resp_matrix[hash_key])
        for index, log in enumerate(self.resp_matrix[hash_key]):
            self.resp_matrix[hash_key][index] /= shift_sum


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
    for pos in range(pwm.length):
        col_count = {nucl: 0 for nucl in pwm.alphabet}
        for seq in seqs:
            hash_key = tools.get_hash(seq)
            for offset in range(len(seq) - pwm.length + 1):
                snippet = tools.snip(seq, pwm.length, offset)
                col_count[snippet[pos]] += rpm.resp_matrix[hash_key][offset]
        col_total = sum(col_count.values())
        for nucl in pwm.alphabet:
            col_val = col_count[nucl]
            up_val = (col_val + pwm.beta * p0[nucl]) / (col_total + pwm.beta)
            pwm.update(up_val, nucl, pos)

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
    p0_total = sum(new_p0.values())
    for nucl in pwm.alphabet:
        new_p0[nucl] /= p0_total
    p0.clear()
    p0.update(new_p0)