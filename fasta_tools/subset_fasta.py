#!/usr/bin/env python3

import sys
import pysam

i = sys.argv[2]

with pysam.FastxFile(sys.argv[1]) as fh:
    for entry in fh:
        if i in entry.name:
            seq = ">{}\n{}\n".format(entry.name, entry.sequence)
            print(seq)
            
