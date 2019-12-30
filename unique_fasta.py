#!/usr/bin/env python3

import pysam
import sys
import os


def main(argv):
        if len(argv) != 2:
            sys.stderr.write("Usage: %s <fastafile> <fastaoutfile>\n" % os.path.basename(sys.argv[0]))
            sys.exit(1)
        fasta = sys.argv[1]

        fh = pysam.FastxFile(fasta)        
        fa_out     = open(sys.argv[2], 'w')
        
        unique_seqs = []
        for entry in fh:
            if entry.name in unique_seqs:
                continue
            else:
                seq            = entry.sequence
                out_seq = '>{0}\n{1}\n'.format(entry.name, seq)
                fa_out.write(out_seq)
                unique_seqs.append(entry.name)
        fa_out.close()
        fh.close()
            
            
            
if __name__ == '__main__':
    main(sys.argv[1:])
