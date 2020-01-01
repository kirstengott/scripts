#!/usr/bin/env python3

import pysam
import sys
import os




def rename_fasta_seqs(fasta):
    fh = pysam.FastxFile(fasta)
    fa_out     = open(sys.argv[2], 'w')
    fasta_base = os.path.basename(fasta)
    fa_out_dir = os.path.dirname(sys.argv[2])
    id_map     = fasta_base.split(".")[:-1]
    if len(id_map) > 1:
        id_map = ".".join(id_map)
    else:
        id_map = id_map[0]
    id_map = fa_out_dir + "/" + id_map + "_id_map.txt"
    id_map_out = open(id_map, 'w')
    count=0
    for entry in fh:
        seq            = entry.sequence
        new_id = 'seq' + str(count)
        out_seq = '>{0}\n{1}\n'.format(new_id, seq)
        out_id_map = "{0}\t{1}\n".format(new_id, entry.name)
        fa_out.write(out_seq)
        id_map_out.write(out_id_map)
        count += 1
    id_map_out.close()
    fa_out.close()
    fh.close()

    

def main(argv):
    if len(argv) != 2:
        sys.stderr.write("Usage: %s <fastafile> <fastaoutfile>\n" % os.path.basename(sys.argv[0]))
        sys.exit(1)
    fasta = sys.argv[1]
    rename_fasta_seqs(fasta)



    
if __name__ == '__main__':
    main(sys.argv[1:])    
    
