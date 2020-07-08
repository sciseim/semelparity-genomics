# ################################################################################
# code from https://www.reddit.com/r/bioinformatics/comments/1u8yc7/looking_for_a_script_that_will_split_dna/
# ################################################################################


'''
Splits all sequences within a multi-fasta file into chunks of a specified size.

Fasta header information is retained with each split sequence - its position in
the original is appended to the id. Single-line and multi-line fasta files are
supported. Prints to stout, so pipe to a file to store the result.

Usage:
python splitter.py <filename> <chunksize>
python splitter.py myfile.fa 100
'''

from __future__ import print_function
from sys import argv, version

if version[0] == '2':
    from itertools import izip_longest as zl
else:
    from itertools import zip_longest as zl

chunksize = int(argv[2])

def writeseq(header, seq):
    for i, chunk in enumerate(zl(*[iter(seq)]*chunksize, fillvalue='')):
        print(header + '_{}bp'.format(i*chunksize))
        print(''.join(chunk))

with open(argv[1]) as f:
    header, seq = f.readline().rstrip(), ''
    for l in f:
        if l[0] != '>':
            seq += l.rstrip()
        else:
            writeseq(header, seq)
            header, seq = l.rstrip(), ''
    writeseq(header, seq)
