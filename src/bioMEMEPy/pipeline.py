def meme(fasta, alphabet, motif_l, loops):
    # Errors
    if not isinstance(motif_l, int):
        raise TypeError('Length of the motif must be int.')
    if not isinstance(loops, (int or bool)):
        raise TypeError('Amount of loops must be either int or None.')