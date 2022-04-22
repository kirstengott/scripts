#!/usr/bin/env python3

import sys, os
from Bio import SeqIO
from Bio.SeqIO import AbiIO


def main(argv):
    if len(argv) < 2:
        sys.stderr.write('''Usage: %s <abi-forward> <abi-reverse>
        This script takes the forward and reverse abi files from sanger sequencing
        trims them, and joins them.
        <abi-forward>   the forward sequence file
        <abi-reverse>   the reverse sequence file \n''' % os.path.basename(sys.argv[0]))
        sys.exit(1)

    forw = SeqIO.parse(sys.argv[1], "abi-trim")
    rev  = SeqIO.parse(sys.argv[2], "abi-trim")

    name = os.path.splitext(os.path.basename(sys.argv[1]))[0] + ':' + os.path.splitext(os.path.basename(sys.argv[2]))[0]
        

    seq1 = ''
    seq2 = ''
    for seq in forw:
        seq1 = seq.seq
    for seq in rev:
        seq2 = seq.seq.reverse_complement()


    print(">{}\n{}{}".format(name, seq1, seq2))


    
if __name__ == '__main__':
    main(sys.argv[1:])
