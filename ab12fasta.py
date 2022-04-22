#!/usr/bin/env python3

import sys, os
from Bio import SeqIO
from Bio.SeqIO import AbiIO

forw = SeqIO.parse(sys.argv[1], "abi-trim")
rev  = SeqIO.parse(sys.argv[2], "abi-trim")

name = os.path.splitext(os.path.basename(sys.argv[1]))[0] + ':' + os.path.splitext(os.path.basename(sys.argv[2]))[0]


seq1 = ''
seq2 = ''
for seq in forw:
    seq1 = seq.seq
for seq in rev:
    seq2 = seq.seq


print(">{}\n{}{}".format(name, seq1, seq2))

