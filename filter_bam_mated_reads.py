#!/usr/bin/env python

import pysam
import sys
import os

## pull in arguments
file_in = sys.argv[1]
out_dir = os.path.dirname(file_in)

## define output files
out_f_mated   = out_dir + "/allpaired.bam"
#out_f_mated2   = out_dir + "/allpaired_R2.fastq.gz"
out_f_unmated  = out_dir + "/allunmated.bam"

## open input file
samfile             = pysam.AlignmentFile(file_in, "rb")

## open output files
mated_aligned_read = pysam.AlignmentFile(out_f_mated, "wb", template = samfile)
unmated_aligned_read = pysam.AlignmentFile(out_f_unmated, "wb", template = samfile)

#format_fastq = "@{0}\n{1}\n+\n{2}\n".format

for read in samfile.fetch():
         if read.is_proper_pair:
                  mated_aligned_read.write(read)
         elif read.mate_is_unmapped:
                  unmated_aligned_read.write(read)
                  # unmated_aligned_read.write(
                  #          format_fastq(read.query_name,
                  #                       read.query_sequence,
                  #                       ''.join(
                  #                                map(lambda x: chr(x+33), read.query_qualities))))




mated_aligned_read.close()
#mated_aligned_read2.close()
samfile.close()
unmated_aligned_read.close()


