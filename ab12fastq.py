#!/usr/bin/env python3

import sys
from Bio import SeqIO
from Bio.SeqIO import AbiIO

forw = SeqIO.parse(sys.argv[1], "abi-trim")
rev  = SeqIO.parse(sys.argv[2], "abi-trim")



seq1 = ''
seq2 = ''
for seq in forw:
    seq1 = seq.seq
for seq in rev:
    seq2 = seq.seq


print("{}{}".format(seq1, seq2))

