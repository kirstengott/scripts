#!/usr/bin/env python3

import pysam
import os
import sys


input_dir = sys.argv[1]
input_fa_dir = sys.argv[2]
busco_out_dir = sys.argv[3]
os.mkdir(busco_out_dir)



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



dirs = os.listdir(input_dir)


all_busco2gene = {}
key_counts = {}
total_db_length = 1

for d in dirs:
    if 'run_' in d:
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
    
buscos_pass = []

## determine which buscos are in every genome
for i in key_counts:
    if key_counts[i] == total_db_length-1:
        buscos_pass.append(i)
#print(total_db_length)


#exit()


buscos_pass_sequences = {}

for i in all_busco2gene.keys():
    busco_dictionary = all_busco2gene[i]
    fasta_reference = input_fa_dir + "/" + i
    keep_seqs = [busco_dictionary[x] for x in buscos_pass]
    sequences = fetch_fasta(fasta_reference, keep_seqs)
    for b, y in zip(buscos_pass, sequences):
        out_file = busco_out_dir + "/" + b + ".fa"
        out_f = open(out_file, 'a')
        out_f.write(">" + i + "\n") ## write out the sequence with the original FASTA file as the id (for tree building) 
        out_f.write(y + "\n")

        


## iterate over the fasta and pull out the ids we want


