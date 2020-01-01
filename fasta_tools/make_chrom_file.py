#!/usr/bin/python

from Bio import SeqIO
import argparse


parser = argparse.ArgumentParser(description="")
parser.add_argument("-i", help = "path to input FASTA")
parser.add_argument("-o", help = "output path")
args = parser.parse_args()


in_fasta = open(args.i)

records = SeqIO.parse(in_fasta, 'fasta')

out_file = open(args.o, 'w')

for r in records:
    out_file.write(r.id + '\t' + str(len(r.seq)) + '\n')

in_fasta.close()
out_file.close()
