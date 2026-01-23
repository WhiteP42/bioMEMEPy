import logging
import time
from . import oops
from . import zoops
from . import tools
from .tools import BasePWM as PWM


# Logger call:
logger = logging.getLogger(__name__)


def meme(fasta, alphabet, model, m_length, motif_num=1, top_val=0.5, seed_limit=5000, threshold=1e-5, max_iter=200,
         emp_background=False, lazy=False):
    start_time = time.perf_counter()
    if not isinstance(m_length, int):
        raise TypeError('Length of the motif must be int.')
    seqs = tools.extract(fasta)

    #Call mode
    _models = {'oops': oops,
               'zoops': zoops} # ZOOPS and ANR models to be added here.

    #Validate model call:
    if model not in _models:
        raise ValueError(f'Model {model} not supported.')
    else:
        model = _models[model]

    #Gather seeding candidates (universal):
    logger.debug('Gathering seed candidates.')
    if tools.snip_count(seqs, m_length) > seed_limit:
        seed_seqs = tools.gather(seqs, m_length, seed_limit)
    else:
        seed_seqs = tools.gather(seqs, m_length)
    logger.debug('Done!')
    logger.debug(f'Total seeds: {len(seed_seqs)}.')

    #Seeding process:
    top_candidate = None
    logger.debug(f'Beginning seeding process.')
    for index, snip in enumerate(seed_seqs):
        logger.debug(f'Testing {snip} ({index + 1}/{len(seed_seqs)}).')
        if emp_background:
            p0 = tools.p0_gen(seqs, alphabet)
        else:
            p0 = {nucl: 0.25 for nucl in alphabet}
        pwm = PWM(snip, alphabet, m_length, top_val)
        rpm = model.RPM(m_length) # <-- RPM is model dependant.

        # Seed according to seeding mode:
        model.e_step(pwm, rpm, seqs, p0)
        if lazy:
            raise NotImplementedError # Later implementation
        else:
            model.m_step(pwm, rpm, seqs, p0)
            rpm = model.RPM(m_length)
            model.e_step(pwm, rpm, seqs, p0)

        log_like = rpm.log_like
        if top_candidate is None or log_like > top_candidate[1]:
            logger.debug(f'Log-like computed: {log_like}. Replacing top candidate.')
            top_candidate = [snip, log_like]
        logger.debug(f'Log-like computed: {log_like}. Top candidate is {top_candidate[1]}.')

    # Run EM until convergence:
    logger.debug('Beginning EM...')
    seed = top_candidate[0]
    if emp_background:
        p0 = tools.p0_gen(seqs, alphabet)
    else:
        p0 = {nucl: 0.25 for nucl in alphabet}
    pwm = PWM(seed, alphabet, m_length, top_val)
    converged = False
    prev_loglike = None
    iterations = 0
    while not converged:
        iterations += 1
        logger.debug(f'Iteration {iterations}.')
        rpm = model.RPM(m_length)  # <-- RPM has to be cleared every loop.
        model.e_step(pwm, rpm, seqs, p0)
        model.m_step(pwm, rpm, seqs, p0)
        log_like = rpm.log_like
        logger.debug(f'Log-like is {log_like}.')
        if prev_loglike is None or (log_like - prev_loglike) > threshold:
            prev_loglike = log_like
        elif (log_like - prev_loglike) <= threshold:
            converged = True
        if iterations >= max_iter:
            converged = True

    # Generate consensus from PWM:
    consens = tools.consensus(pwm.matrix, m_length, alphabet)

    total_time = time.perf_counter() - start_time
    logger.info(f'Runtime: %.3f s.', total_time)

    #Return PWM
    return pwm.matrix, consens