#!/usr/bin/env python3

import pysam
import sys
import re
import os
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

def main(argv):
    if len(argv) < 2:
        sys.stderr.write('''Usage: %s <fastafile> <gff> <cds_from_genomic_out>
        <fastafile> : the fasta file to extract cds region from
        <gff> : gff file with cds regions specified
        <cds_from_genomic_out> : output file\n''' % os.path.basename(sys.argv[0]))
        sys.exit(1)

    fasta = sys.argv[1]
    gff   = sys.argv[2]
    outfile = sys.argv[3]
    gff_cds = {}

    with open(gff, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            line = line.strip().split()
            if line[2] in ['cds', 'CDS']:
                note = line[8].split(";")
                if 'LGSR00000000' in fasta:
                    parent_ind = [note.index(x) for x in note if 'Name' in x][0]
                    parent = re.sub("Name=", "", note[parent_ind])
                else:
                    parent_ind = [note.index(x) for x in note if 'Parent' in x][0]
                    parent = re.sub("Parent=", "", note[parent_ind])
                id_ind = [note.index(x) for x in note if 'ID=' in x][0]
                name = re.sub("ID=", "", note[id_ind])
                chrom = line[0]
                start = int(line[3])
                stop = int(line[4]) 
                strand = line[6]
                region = '{}:{}-{}'.format(chrom, start, stop)
                if parent not in gff_cds.keys():
                    gff_cds[parent] = [[region, strand]]
                else:
                    gff_cds[parent].append([region, strand])


    fasta_in = pysam.FastaFile(fasta)
    out_f = open(outfile, 'w')
    for mrna in gff_cds:
        all_cds = gff_cds[mrna]
        full_cds = ''
        for cds in all_cds:
            strand  = cds[1]
            seq     = fasta_in.fetch(region = cds[0])
            if strand == "-":
                seq = Seq(seq, generic_dna)
                seq = seq.reverse_complement()
            full_cds += seq
        outstring = ">{}\n{}\n".format(mrna, full_cds)
        out_f.write(outstring)
    out_f.close()
    fasta_in.close()
    




if __name__ == '__main__':
    main(sys.argv[1:])
