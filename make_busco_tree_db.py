#!/usr/bin/env python3

import pysam
import os
import sys
import re

def get_busco_completes(busco_file):
    ## define a dictionary where the key is the BUSCO id and the value is the sequence name
    busco_completes = {}
    for line in busco_file:
        if line.startswith("#"):
            pass
        else:
            entries = line.strip().split("\t")
            if entries[1] == 'Complete':
                busco_completes[entries[0]] = entries[2]
    return(busco_completes)

def fetch_fasta(fasta_reference, sequences):
    ## pull out individual fasta records
    fh = pysam.FastaFile(fasta_reference)
    seqs_all = []
    for rec in sequences:
        seq = fh.fetch(reference = rec)
        seqs_all.append(seq)
    return(seqs_all)
    fh.close()

def write_treedb_diction2fasta(matches_dictionary, fasta_dir,  output_directory, keep_sequences = False):
    if not os.path.exists(output_directory):
        os.mkdir(output_directory)
        print('Made directory:', output_directory)
    if os.path.exists('seqid_map.txt'):
        os.remove('seqid_map.txt')
    id_map = open('seqid_map.txt', 'a')
    ind = 1
    for i in matches_dictionary.keys():
        sub_dictionary = matches_dictionary[i]
        fasta_reference = fasta_dir + "/" + i
        if keep_sequences:
            keep_seqs = [sub_dictionary[x] for x in keep_sequences]
            names_scheme = keep_sequences
        else:
            keep_seqs = sub_dictionary.values()
            names_scheme = sub_dictionary.keys()
        ## fetch the kept sequences
        sequences = fetch_fasta(fasta_reference, keep_seqs)
        for b, y in zip(names_scheme, sequences):
            out_file = output_directory + "/" + re.sub('.*\|', '', b) + ".fa"
            out_f = open(out_file, 'a')
            out_f.write(">seq" + str(ind) + "\n") ## write out the sequence with the original FASTA file as the id (for tree building) 
            out_f.write(y + "\n")
            out_f.close()
        id_map.write(i + "\tseq" + str(ind) + "\n")
        ind += 1
    ## finish out the function    
    id_map.close()



def get_shared_ids(key_counts, total_db_length):
    id_pass = []
    for i in key_counts:
        if key_counts[i] == total_db_length:
            id_pass.append(i)
    return(id_pass)



def main(argv):
    if len(argv) != 3:
        sys.stderr.write("Usage: %s <in.busco.output.dir> <in.reference.fasta.dir> <output.dir>\n" % os.path.basename(sys.argv[0]))
        sys.exit(1)
    input_dir     = sys.argv[1]
    input_fa_dir  = sys.argv[2]
    busco_out_dir = sys.argv[3]

    ## read in BUSCO data and create a dictionary of all of the results
    dirs = os.listdir(input_dir)
    all_busco2gene  = {}
    key_counts      = {}
    total_db_length = 1
    for d in dirs:
        if 'run_' in d:
            print(d)
            total_db_length += 1
            f = d[4:] ## subset the file name to just the GCA accession to reference the source fasta file with
            b_file = input_dir + "/" + d +"/full_table_" + f + ".tsv"
            summary_file = input_dir + "/" + d +"/short_summary_" + f + ".txt"
            try:
                sum_file = open(summary_file, 'r')
                lines = [x.strip() for x in sum_file]
                if not '3.1.0' in lines[0]:
                    print('BUSCO version not 3.1.0, program may fail')
            except:
                print('No required records in ' + input_dir + "/" + d  + " removing")
                os.rmdir(input_dir + "/" + d )
            try:
                busco_file = open(b_file, 'r')
                b = get_busco_completes(busco_file)
                all_busco2gene[f] = b ## save the dictionary and source name
                if not key_counts:
                    for i in b.keys():
                        key_counts[i] = 1
                else:
                    for i in b.keys():
                        if i in key_counts:
                            key_counts[i] += 1
                        else:
                            pass
            except:
                print('No required records in ' + input_dir + "/" + d  + " removing")
                os.rmdir(input_dir + "/" + d )
        else:
            pass
    

    ## determine which buscos are in every genome
    buscos_pass = get_shared_ids(key_counts = key_counts, total_db_length = total_db_length-1)

    ## write out the sequences for our tree
    write_treedb_diction2fasta(matches_dictionary = all_busco2gene,
                               fasta_dir = input_fa_dir,
                               keep_sequences = buscos_pass,
                               output_directory = busco_out_dir)


if __name__ == '__main__':
    main(sys.argv[1:])
