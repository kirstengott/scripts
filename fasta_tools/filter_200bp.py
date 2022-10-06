import pysam
import sys
import os


fasta = sys.argv[1]


fh = pysam.FastxFile(fasta)

for entry in fh:
    if len(entry.sequence) < 200:
        pass
    else:
        out_seq = '>{0}\n{1}\n'.format(entry.name, entry.sequence)
        print(out_seq)
    

