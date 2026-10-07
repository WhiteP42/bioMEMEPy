"""
Currently OOPS generation only, for testing purposes
"""
import random

motif = 'ATAACGTATTACGGCGGCAT'
alphabet = ['A', 'T', 'C', 'G']
seq_len = 500
n_seq = 3000
break_range = (3, 10)

for seq in range(n_seq):
    input_seq = []
    motif_inj = random.randint(0, seq_len - len(motif))
    sec_start = seq_len - (motif_inj + len(motif))
    for pos in range(motif_inj):
        input_seq.append(random.choice(alphabet))
    input_seq.append(motif)
    for pos in range(sec_start):
        input_seq.append(random.choice(alphabet))
    input_seq = ''.join(input_seq)
    with open('oopstest.fasta', 'a') as f:
        f.write(f'>seq{seq}: m_pos = {motif_inj}, motif = {motif}\n'
                f'{input_seq}\n')