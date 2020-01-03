#!/usr/bin/env python3

import pysam
import sys
import os




def rename_fasta_seqs(fasta, fasta_out):
    '''Renames of fasta file with 'seqN' where N is a sequential number
       and writes out a file with extenstion '_id_map.txt' with mappings
       to the original names fasta sequence names.
       fasta: path to the fasta file
       fasta_out: path to the output fasta file
    '''
    fh = pysam.FastxFile(fasta)
    fa_out     = open(fasta_out, 'w')
    fasta_base = os.path.basename(fasta)
    fa_out_dir = os.path.dirname(fasta_out)
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
    fasta_out = sys.argv[2]
    rename_fasta_seqs(fasta = fasta, fasta_out = fasta_out)




    
if __name__ == '__main__':
    main(sys.argv[1:])    
    
