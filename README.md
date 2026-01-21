# bioMEMEPy (v0.1.0)

>**WARNING: THIS PACKAGE IS STILL IN EARLY ALPHA! BUGS ARE EXPECTED.**

`bioMEMEPy` will be pip-installable Python package for motif discovery using **OOPS**, **ZOOPS** and **ANR**.

## Current status of the project:

- OOPS mode: **Available**, pending optimization.
  - Algorithm fully implemented.
  - Fasta file generator available in `tests/fasta/oops/` for algorithm testing along with `crp0.fasta`.
  - Only allows for one motif to be discovered; multiple discovery will be added later.
  - Doesn't use ambiguity codes yet, planned feature.
- ZOOPS mode: *Planned*
- ANR mode: *Planned*

## How to use bioMEMEPy:

Install the package by using `pip install git+https://github.com/WhiteP42/bioMEMEPy.git`.

The main command line is `meme(fasta, alphabet, model, m_length)`, all explained below:<br>
`fasta`: Location of the fasta file.<br>
`alphabet`: Nucleotides or aminoacids your fasta file sequences use.<br>
`model`: The model you plan on using. Currently, only OOPS is available.<br>
`m_length`: The motif length. Because ambiguity codes are not yet implemented, it is recommended to keep this value as
close as possible to the actual motif length.

The command will return `pwm`, a dictionary with the frequencies for all positions of the motif; and `consensus`, a
string containing the consensus motif. The motif will be returned with absolute nucleotides; this means no ambiguity
will be displayed. You can, however, use the PWM provided to analize it. Ambiguity codes are a feature planned for
later.

You also have, if required, a nucleotide or aminoacid list, which you can access by calling `bioMEMEPy.mnm.dna`, 
`bioMEMEPy.mnm.rna` or `bioMEMEPy.mnm.aa`.

There are also optional parameters, with default values assigned to, you can configure if you wish:<br>
`motif_num=1`: This parameter is not implemented yet; it will allow you to set the number of motifs MEME should find.<br>
`top_val=0.5`: Sets the value that the present nucleotide in a newly generated PWM will have per position.<br>
`seed_limit=5000`: Sets the number of motifs that will be evaluated during the seeding process. Increase for rare motif,
decrease for better performance. **Never decrease below 1000!**<br>
`threshold=1e-5`: Value which sets convergence requirements for the EM algorithm.<br>
`max_iter`: Maximum number of EM iterations.<br>
`d_background=False`: This mode permits the generation of dynamic starting bacground proportions, which counts all
characters in the dataset. Leave as `False` for DNA and RNA; switch to `True` for proteins.<br>
`lazy=False`: Optimization trick that disables M-step when evaluating seeds, which makes it so it only evaluates the
log probabilities of the starting sequence and not the one derived from one full EM loop, lowering seed quality.
Switch to `True` for a huge performance boost at the cost of accuracy. *Recomended for large datasets.*

Example (basic):<br>
`pwm, consensus = bioMEMEPy.meme('fastas/crp0.fasta', bioMEMEPy.mnm.dna, 'oops', 16)`

Example (changes optional parameter `lazy` and `seed_limit`):<br>
`pwm, consensus = bioMEMEPy.meme('fastas/crp0.fasta', bioMEMEPy.mnm.dna, 'oops', 16, seed_limit=10000, lazy=True)`