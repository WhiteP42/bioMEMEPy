from . import handler

# OOPS model
def oops(fasta, alphabet, motif_l, loops):
    # Errors
    if not isinstance(motif_l, int):
        raise TypeError('Length of the motif must be int.')
    if not isinstance(loops, (int or bool)):
        raise TypeError('Amount of loops must be either int or None.')

    seqs = handler.extract(fasta)
    pwm = handler.PWM(alphabet, motif_l)
    # Initialization
    pwm.init()
    finish = False
    loop_count = 0
    while loop_count != loops or not finish:
        loop_count += 1
        raise NotImplementedError