#!/usr/bin/env python3

import re
import sys
import pysam
import os
import math






def main(argv):
    if len(argv) != 2:
        sys.stderr.write("Usage: %s <fastafile> <n_seq_file>\n" % os.path.basename(sys.argv[0]))
        sys.exit(1)
    file_n = sys.argv[1]
    fasta_file = pysam.FastxFile(file_n)
    max_num_seqs = int(sys.argv[2])


    length = 0

    for x in fasta_file:
        length += 1

    n_files = math.ceil(length / max_num_seqs)


    count = 1
    file_init = 1

    f_out = "sequences_" + str(file_init) + ".fa"

    f_o = open(f_out, 'w')

    fasta_file = pysam.FastxFile(file_n)



    for entry in fasta_file:
        seq            = entry.sequence
        out = '>{0}\n{1}\n'.format(entry.name, seq)
        if count <= max_num_seqs:
            f_o.write(out)
            count += 1
        else:
            count = 1
            file_init += 1
            f_o.close()
            f_out = "sequences_" + str(file_init) + ".fa"
            f_o = open(f_out, 'w')
            f_o.write(out)
        
if __name__ == '__main__':
    main(sys.argv[1:])
