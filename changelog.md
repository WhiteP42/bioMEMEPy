***v0.1.3***

**CORE CHANGES**
- Loosened Python requirements to **v3.9** to better fit the scientific community.
- Fixed the selection algorithm so it doesn't stall during the collection of k-mers when ```amount > 0```, making the
algorithm faster.

**TEST CHANGES**
- Changed ```fasta_gen.py``` to be able to generate files making them miss a portion of the motif dinamically. This will
enhance the amount of tests we can perform by generating mock sequence data befor jumping to real-life data, which we
will do when the ZOOPS algorithm is implemented.

*Changelogs begun being kept at version 0.1.3*