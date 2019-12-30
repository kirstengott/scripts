#!/usr/bin/env python

import pysam
import sys
import os
import re
import subprocess


def write_read(fastq, read):
             """
    Write read to open FASTQ file.
    """
             info = {'index': int(not read.is_read1) + 1,
                     'name':  read.qname}
             if read.is_reverse:
                      info.update({'quality':  read.qual[::-1],
                                   'sequence': reverse_complement(read.seq)})
             else:
                      info.update({'quality':  read.qual,
                                   'sequence': read.seq})
                      fastq.write('@{name}/{index}\n{sequence}\n+\n{quality}\n'.format(**info))

def main(argv):
         ## pull in the bam file
         if len(argv) != 4:
                  sys.stderr.write("Usage: %s <bam_file> <contaminant_reference_fa> <output_directory> <s>\n" % os.path.basename(sys.argv[0]))
                  sys.exit(1)
                  
         bam_file = sys.argv[1]
         reference = sys.argv[2] ## the name of the file with the reference contaminant that was aligned to
         ## open input file
         bamfile             = pysam.AlignmentFile(bam_file, "rb")

         ref_fasta     = pysam.FastaFile(reference)
         reference_ids = ref_fasta.references

         ## define output files
         out_dir = sys.argv[3]


         if os.path.abspath(out_dir) == os.path.abspath(os.getcwd()):
             out_f_mated = "filtered_" + os.path.basename(bam_file)
             out_base = re.sub('.bam', '', out_f_mated)          ## output fastqfile
             contaminant_reads = "contaminant_" + os.path.basename(bam_file)
         else:
             out_f_mated   = out_dir + "/" + os.path.basename(bam_file)          ## output bam file
             out_base = re.sub('.bam', '', os.path.basename(bam_file))          ## output fastqfile
             contaminant_reads = "contaminant_" + os.path.basename(bam_file)
         ## open output files
         contam_reads = pysam.AlignmentFile(contaminant_reads, 'wb', template = bamfile)
         mated_unaligned_read = pysam.AlignmentFile(out_f_mated, "wb", template = bamfile)
         # fastq_r1_o = open(fastq_r1, 'w')
         # fastq_r2_o = open(fastq_r2, 'w')
         

         
         # read1 = None
         # read2 = None
         
         for read in bamfile:
                      ## logic to test if the read or its paired end map to a contaminant reference id
             if sys.argv[4] == 's':
                          if read.reference_name in reference_ids:
                                       contam_reads.write(read)
                                       continue ## if either read maps to contaminant, continue to the next iteration
                          else:
                                       mated_unaligned_read.write(read) ## otherwise write the pairs out

                          

             else:
                          if read.reference_name in reference_ids or read.next_reference_name in reference_ids:
                                       contam_reads.write(read)
                                       continue ## if either read maps to contaminant, continue to the next iteration
                          else:
                                       mated_unaligned_read.write(read) ## otherwise write the pairs out
         mated_unaligned_read.close()
         bamfile.close()
         contam_reads.close()

         if sys.argv[4] == 's':
                      fastq_out = out_dir + "/" + out_base + ".fastq"
                      subprocess.run(["bedtools", "bamtofastq", "-i", out_f_mated, "-fq", fastq_out,])
                      subprocess.run(["gzip", fastq_out])

         else:
                      fastq_r1 = out_dir + "/" + out_base + "_R1.fastq"
                      fastq_r2 = out_dir + "/" + out_base + "_R2.fastq"
                      subprocess.run(["bedtools", "bamtofastq", "-i", out_f_mated, "-fq", fastq_r1, '-fq2', fastq_r2])
                      subprocess.run(["gzip", fastq_r1, fastq_r2])

         os.remove(out_f_mated)


if __name__ == '__main__':
    main(sys.argv[1:])
