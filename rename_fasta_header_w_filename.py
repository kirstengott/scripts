#!/usr/bin/env python3

import re
import sys
import pysam


def main(argv):
    if len(argv) != 1:
        sys.stderr.write("Usage: %s <in.fastafile>\n" % os.path.basename(sys.argv[0]))
        sys.exit(1)
    file_n = sys.argv[1]
    fasta_file = pysam.FastxFile(file_n)
    for entry in fasta_file:
        seq            = entry.sequence
        out = '>{0};{1}\n{2}'.format(entry.name, file_n, seq)
        print(out)

    
if __name__ == '__main__':
    main(sys.argv[1:])
