1. Make an artificial DNA sequence of 200-400b in length, in which to simulate 3-4 transposable elements.
2. Implement a software application to detect the positions of these transposable elements (start, end) within the created DNA sequence.

Download a FASTA file from NCBI that contains an entire genome (bacterial / plant genome up to 6mb (mega-bases??)).
Change your application such that it may find new transposoms across the entire genome. 

Note: The inverted repeats this time are unknown.
Proposal: One way to solve this task will be for you to use 2 different sliding windows S1 and S2. 
->The first sliding window (S1) should slide across the entire genome with a step of 1, but each new position of S1, the S2 sliding window should search for the second pattern up to the end of the genome. Note that the second pattern should be the content found in S1, inverted, and what you have to search is the complementary part of this inversion.

Note: The minimum length between S1 and S2 must be at least 100 bases. (S2 starts 100 bases after S1 and starts until the end)