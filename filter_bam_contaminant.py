#!/usr/bin/env python

import pysam
import sys
import os

## pull in arguments
file_in = sys.argv[1]
out_dir = os.path.dirname(file_in)

## define output files
out_f_mated   = out_dir + "/contaminant_filtered.bam"


## open input file
samfile             = pysam.AlignmentFile(file_in, "rb")
## open output files
mated_aligned_read = pysam.AlignmentFile(out_f_mated, "wb", template = samfile)

ref_fasta = pysam.FastaFile(sys.argv[2])

reference_ids = ref_fasta.references


for read in samfile:
         if read.reference_name in reference_ids or read.next_reference_name in reference_ids:
                  pass
         else:
                  mated_aligned_read.write(read)



mated_aligned_read.close()
samfile.close()


