from . import handler

#OOPS model

def oops(fasta, loops, alphabet):
    seqs = handler.extract(fasta)
    pwm = handler.PWM(alphabet, len(seqs))
    # Initialization
    pwm.init()