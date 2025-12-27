def meme(fasta, alphabet, motif_l, loops, model):
    # Errors
    if not isinstance(motif_l, int):
        raise TypeError('Length of the motif must be int.')
    if not isinstance(loops, int):
        raise TypeError('Amount of loops must be either int or None.')
    if model not in ['oops', 'zoops', 'anr']:
        raise ModuleNotFoundError('Mode not supported.')