Take an arbitrary DNA sequence from NCBI (bacterial genome). 
Use 5 restriction enzymes (name, recognized site, cleavage site)

A restriction enzyme is a protein that does chemical work inside the cell.
It is a special type of enzyme, it is called restriction because it cuts something (cleavages)
The restriction enzymes are specific for DNA sequences

Wherever they see the pattern, they link to it and they start to cut the DNA in a specific position inside the pattern, or outside the pattern.


EcoRI   5'GAATTC    5'---G      ^   AATTC---3'
        3'CTTAAG    3'---CTTAA  ^   G---5'

BamHI   5'GGATCC    5'---G      ^   GATCC---3'
        3'CCTAGG    3'---CCTAG  ^   G---5'
HindIII 5'AAGCTT    5'---A      ^   AGCTT---3'
        3'TTCGAA    3'---TTCGA  ^   A---5'
TagI    5'TCGA      5'---T      ^   CGA---3'
        3'AGCT      3'---AGC    ^   T---5'
HaeIII  5'GGCC      5'---GG     ^   CC---3'
        3'CCGG      3'---CC     ^   GG---5'


Input of the implementation:
1. The recognized sequence for each restriction enzyme
2. The DNA sequence to be digested

Output of the implementation:
1. Number of cleavages for each restriction enzyme
2. Cleavage position and length of fragments
3. A simulation of the electrophoresis gel based on the number of restriction enzymes used