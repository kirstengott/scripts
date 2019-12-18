#!/usr/bin/env python3

import re
import sys
import pysam


file_n = sys.argv[1]
fasta_file = pysam.FastxFile(file_n)
for entry in fasta_file:
    seq            = entry.sequence
    out = '>{0}\n{1}'.format(entry.name, seq)
    print(out)
