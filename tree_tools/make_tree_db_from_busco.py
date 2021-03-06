#!/usr/bin/env python3

import pysam
import os
import sys, shutil
import re
from pathlib import Path




def get_busco_completes(busco_file):
    '''
    Define a dictionary where the key is the BUSCO id and the value is the sequence name
    Only include BUSCO ids that are 'Complete'
    busco_file is the 'full_output' busco output file
    '''
    busco_completes = {}
    dupes = set()
    for line in busco_file:
        if line.startswith("#"):
            pass
        else:
            entries = line.strip().split("\t")
            if entries[1] == 'Complete':
                if entries[0] in busco_completes:
                    dupes.add(entries[0])
                busco_completes[entries[0]] = entries[2]

    for x in dupes:
        busco_completes.pop(x)## remove duplicate buscos (new for v4)
    return(busco_completes)

def fetch_fasta(fasta_reference, sequences):
    '''
    Pull out individual fasta records
    '''
    fh = pysam.FastaFile(fasta_reference)
    seqs_all = []
    for rec in sequences:
        seq = fh.fetch(reference = rec)
        seqs_all.append(seq)
    return(seqs_all)
    fh.close()

def write_treedb_diction2fasta(matches_dictionary,
                               fasta_dir,
                               output_directory,
                               keep_sequences = False):
    ''' 
    For all BUSCO ids, write out a FASTA file with proteins sequences pertaining.
    Create a naming scheme for genomes suitable for RaXML with a file encoding 
    the original genomes nameshoused in 'seqid_map.txt'
    matches_dictionary: dictionary containing BUSCO to genome to gene id mappings
    fasta_dir: directory containing fasta files for genes
    output_directory: where to store the resulting fasta files
    keep_sequences: subset of sequences to output
    ''' 
    if not os.path.exists(output_directory):
        os.mkdir(output_directory)
        print('Made directory:', output_directory)
    if os.path.exists('id_map.txt'):
        os.remove('id_map.txt')
    id_map = open('id_map.txt', 'a')
    ind = 1
    fasta_files = [os.path.basename(x) for x in Path(fasta_dir).rglob('*a')]
    for i in matches_dictionary.keys():
        sub_dictionary = matches_dictionary[i]
        fa_match_file = [x for x in fasta_files if i in x][0]
        fasta_reference = os.path.join(fasta_dir, fa_match_file)
        if keep_sequences:
            keep_seqs = [sub_dictionary[x] for x in keep_sequences]
            names_scheme = keep_sequences
        else:
            keep_seqs = sub_dictionary.values()
            names_scheme = sub_dictionary.keys()
        ## fetch the kept sequences
        sequences = fetch_fasta(fasta_reference, keep_seqs)
        ## write out the fasta files
        for b, y in zip(names_scheme, sequences):
            out_file = output_directory + "/" + re.sub('.*\|', '', b) + ".fa"
            out_f = open(out_file, 'a')
            out_f.write(">seq" + str(ind) + "\n") ## write out the sequence with the original FASTA file as the id (for tree building) 
            out_f.write(y + "\n")
            out_f.close()
        id_map.write('seq' + str(ind) + "\t" + i + "\n")
        ind += 1
    ## finish out the function    
    id_map.close()

    
def get_shared_ids(key_counts, total_db_length):
    '''
    Return the BUSCO ids that are present in all species
    '''
    id_pass = []
    for i in key_counts:
        if key_counts[i] == total_db_length:
            id_pass.append(i)
    return(id_pass)


def busco_3_1_0(b_out_files, input_fa_dir, busco_out_dir):
    ''' 
    This function runs the program method for output from BUSCO version 3.1.0,
    which changed it's output to no longer include fasta files of each busco sequences
       b_out_files: a list of the relative paths to all busco output files included
    '''
    all_busco2gene  = {}
    key_counts      = {}
    total_db_length = 1
    
    for f in b_out_files:
        b_file = b_out_files[f]
        if 'run_' in f:
            f = re.sub("run_", "", f)
        total_db_length += 1
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
    ## determine which buscos are in every genome
    buscos_pass = get_shared_ids(key_counts = key_counts, total_db_length = total_db_length-1)

    ## write out the sequences for our tree
    write_treedb_diction2fasta(matches_dictionary = all_busco2gene,
                               fasta_dir = input_fa_dir,
                               keep_sequences = buscos_pass,
                               output_directory = busco_out_dir)

def find_all(name, path):
    result = []
    for root, dirs, files in os.walk(path):
        if name in files:
            result.append(os.path.join(root, name))
    return result

import glob

def busco_4_0(b_out_files, input_dir, busco_out_dir, seq_type = 'fna'):
    ''' 
    This function runs the program method for output from BUSCO version 4.0,
       b_out_files: a list of the relative paths to all busco output files included
    '''
    all_busco2gene  = {}
    key_counts      = {}
    total_db_length = 0

    genomes_keep = []

    c = set()
    for f in b_out_files:
        b_file = b_out_files[f]
        genome = str(b_file).split("/")[1]
        genomes_keep.append(genome)
        single_copies = os.path.join(os.path.dirname(b_file), "busco_sequences/single_copy_busco_sequences/")
        all_sc = os.listdir(single_copies)
        all_s = [x for x in all_sc if seq_type in x]
        b = set(all_s)
        if c == set():
            c = b
        else:
            c = b.intersection(c)

    
    for busco in c:
        outfile = os.path.join(busco_out_dir, busco)
        all_files = find_all(busco, input_dir)
        for seq_file in all_files:
            genome = seq_file.split("/")[1]
            if genome not in genomes_keep:
                continue
            else:
                command = "cat {} | sed -e 's/>.*$/>{}/' >>{}".format(seq_file, genome, outfile)
                os.system(command)


def busco_all(b_out_files, input_dir, busco_out_dir, seq_type = 'fna'):
    ''' 
    This function runs the program method for output from BUSCO version 4.0,
       b_out_files: a list of the relative paths to all busco output files included
    '''
    all_busco2gene  = {}
    key_counts      = {}
    total_db_length = 0

    genomes_keep = []

    c = set()
    remove = set()
    for f in b_out_files:
        b_file = b_out_files[f]
        genome = str(b_file).split("/")[1]
        genomes_keep.append(genome)
        single_copies = os.path.join(os.path.dirname(b_file), "busco_sequences/single_copy_busco_sequences/")
        all_sc = os.listdir(single_copies)
        multi = os.listdir(os.path.join(os.path.dirname(b_file), "busco_sequences/multi_copy_busco_sequences/"))
        miss = os.listdir(os.path.join(os.path.dirname(b_file), "busco_sequences/fragmented_busco_sequences/"))
        remove.update(multi)
        all_s = [x for x in all_sc if seq_type in x]
        b = set(all_s)
        c.update(b)

    len(c)
    c.difference_update(remove)
    len(c)
    for busco in c:
        outfile = os.path.join(busco_out_dir, busco)
        all_files = find_all(busco, input_dir)
        for seq_file in all_files:
            genome = seq_file.split("/")[1]
            if genome not in genomes_keep:
                continue
            else:
                command = "cat {} | sed -e 's/>.*$/>{}/' >>{}".format(seq_file, genome, outfile)
                os.system(command)



def main(argv):
    if len(argv) < 2:
        sys.stderr.write('''Usage: %s <in.busco.output.dir> <output.dir> [in.reference.fasta.dir]
                         in.busco.output.dir: the directory housing BUSCO results for every species being compared, must be the prefix of the fasta file
                         output.dir: directory to create for storing output
                         in.reference.fasta.dir: fasta file directory used for BUSCO input [optional if not using BUSCOv3.1.0]
                         perc_busco: the percent cutoff of BUSCO completeness to include
                         -all: pull all buscos from select genomes rather than an intersection
                         -aa:        pull out amino acid sequences\n''' % os.path.basename(sys.argv[0]))
        sys.exit(1)
    input_dir     = sys.argv[1]
    busco_out_dir = sys.argv[2]
    comp = float(sys.argv[4])
    version = ''
    if os.path.exists(busco_out_dir):
        print(busco_out_dir, 'exists, removing..')
        shutil.rmtree(busco_out_dir)
    os.mkdir(busco_out_dir)
    if '-aa' in sys.argv:
        seq_type = 'faa'
    else:
        seq_type = 'fna'
    ## read in BUSCO data and create a dictionary of all of the results
    dirs = os.listdir(input_dir)
    b_out_files = {}
    for d in dirs:
        try:
            ## try to find the busco output files and open them
            ## exit if they do not exist for a particlar subdirectory
            d_full = os.path.join(input_dir, d)
            b_file       = [x for x in Path(d_full).rglob('*full_table*')][0]
            summary_file = [x for x in Path(d_full).rglob('*short_summary*')][0]
            sum_file = open(summary_file, 'r')
            lines = [x.strip() for x in sum_file]
            sum_file.close()
            if '3.1.0' in lines[0]:
                version = '3.1.0'
            else:
                version = lines[0].strip().split()[-1]
            for line in lines:
                if line.startswith('C'):
                    complete = float(re.sub("C:", "", line.split("%")[0]))
                else:
                    continue
            if complete >= comp:
                b_out_files[d] = b_file
            else:
                print(d, complete, 'exluded')
        except:
            continue
            #print('No files detected in', d, ',excluding from analysis')
    if version not in ['3.1.0', '4.0.beta1']:
        print('Version:', version, ' not version 3.1.0, may need to add method')
    # if version == '3.1.0':
    #     busco_3_1_0(b_out_files = b_out_files,
    #                 input_fa_dir = sys.argv[3],
    #                 busco_out_dir = busco_out_dir)
    # else:
    if not '-all' in sys.argv:
       busco_4_0(b_out_files = b_out_files,
                  input_dir = input_dir,
                  busco_out_dir = busco_out_dir,
                  seq_type = seq_type)
    else:
        busco_all(b_out_files = b_out_files,
                  input_dir = input_dir,
                  busco_out_dir = busco_out_dir,
                  seq_type = seq_type)



if __name__ == '__main__':
    main(sys.argv[1:])
 
