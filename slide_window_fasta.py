#!/usr/bin/env python

import re
import sys
import pysam

fasta_file = pysam.FastxFile(sys.argv[1])
window     = int(sys.argv[2])
k          = int(sys.argv[3])


def gc_content(sequence):
    gc = 0
    GC = set('GCgc')
    for nt in sequence:
        if nt in GC:
            gc += 1
    return 100.0 * gc / len(sequence)
                                            

def count_kmers(kmer_length, sequence):
    kmer_dict = {}
    for i in range(0, len(sequence) - kmer_length +1):
        # extract the kmer as a substring
        kmer = sequence[i:(i+kmer_length)]
        if kmer in kmer_dict:
            kmer_dict[kmer] += 1
        else:
            kmer_dict[kmer] = 1

    count_distinct = 0
    for k in kmer_dict:
        count_distinct += 1
    return(count_distinct)
                                            

print('contig\twindow\tgc\tdistinct_kmers')

for entry in fasta_file:
    seq            = entry.sequence
    windows        = [seq[i:i+window] for i in range(0, len(seq), window)]
    gc             = [gc_content(x) for x in windows]
    distinct_kmers = [count_kmers(kmer_length = k, sequence = x) for x in windows]
    for i in range(0, len(windows)):
        out = '{0}\t{1}\t{2}\t{3}'.format(entry.name, i, gc[i], distinct_kmers[i])
        print(out)
