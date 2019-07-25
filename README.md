# OneKP
PPRfinder code for analysing 1KP transcriptome data
This version of PPRfinder was developed specifically for transcriptome data and is NOT suitable for genome data  

Input for the program is a) a fasta file of protein sequences and b) domain table output from HMMER 3.2 (http://hmmer.org)

Run hmmsearch from the HMMER 3.2.1 package using the motif profiles developed by Cheng et al. 2016 or Gutmann et al. 2019. The latter are available from doi:10.5061/dryad.XXXX. Recommended hmmsearch settings are ‘--domtblout --noali -E 0.1’. PPRfinder selects the best-scoring non-overlapping chain of motifs from the domain table output of hmmsearch, with a motif score threshold of 0 for all motif types, except SS motifs (score >= 10) and DYW motifs (score >= 30). All motifs on the same strand, whatever frame they were in, are considered in this process. Amino acid sequences from different frames are joined where necessary to connect PPR motif arrays, with ‘X’ residues inserted to maintain the approximate length and to indicate the junction. Protein sequences are filtered such that only sequences with either at least two adjacent PPR motifs and a sum of motif scores >= 40, or a DYW motif, are reported. These are stringent thresholds that we believe eliminate almost all non-PPR sequences.

PPRfinder generates 3 output files:

XXXX.orfs.domt_pprs.fa - fasta format files of PPR protein sequences
XXXX.orfs.domt_beads.txt - tab-separated text file of protein-level info
XXXX.orfs.domt_motifs.txt - tab-separated text file of motif-level info

where the XXXX prefix is derived from the input file names

Columns in beads.txt: protein ID, length in aa, motif arrangement, PPR type (P or PLS), protein score (sum of hmmsearch motif scores)
Columns in motifs.txt: protein ID, start, end, score, sequence, 2nd aa, 5th aa, last aa, motif type
