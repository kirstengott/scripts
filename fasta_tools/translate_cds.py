import sys
import pysam
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna



with pysam.FastxFile(sys.argv[1]) as fh:
    for entry in fh:
        while len(entry.sequence) % 3:
            entry.sequence += 'N'
        seq = Seq(entry.sequence, generic_dna)
        seq = seq.translate()
        outstring = ">{}\n{}\n".format(entry.name, seq)
        print(outstring)

