#!/usr/bin/env python3

import pysam
import sys
import os


def main(argv):
    if len(argv) != 3:
        sys.stderr.write("Usage: %s <fastafile> <fastaoutfile> <id_map_out>\n" % os.path.basename(sys.argv[0]))
        sys.exit(1)
    fasta = sys.argv[1]

    fh = pysam.FastxFile(fasta)



    count=0

    fa_out     = open(sys.argv[2], 'w')
    id_map_out = open(sys.argv[3], 'w')

    for entry in fh:
        seq            = entry.sequence
        new_id = 'seq' + str(count)
        out_seq = '>{0}\n{1}\n'.format(new_id, seq)
        out_id_map = "{0}\t{1}\n".format(new_id, entry.name)
        fa_out.write(out_seq)
        id_map_out.write(out_id_map)
        count += 1

    fa_out.close()
    id_map_out.close()
    fh.close()


    
if __name__ == '__main__':
    main(sys.argv[1:])    
    
