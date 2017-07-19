##################
KMER Trinucleotide
##################

Compare trinucleotide effects for any list of kmers with wieghts.

This algorithm will loop through any file of:

+------+--------+
| kmer | weight |
+------+--------+

For each kmer, it will loop through all available trinucleotides, and then
compare those to **every other kmer at that position**. For each trinucleotide
it will ask how, on *average* does the **central nucleotide** affect the score.
In other words, changing the central nucleotide in any trinucleotide changes
the weight by how much?

This is currently under development, as in its current form it will take several
years to run.
