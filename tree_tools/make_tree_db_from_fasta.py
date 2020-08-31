#!/usr/bin/env python
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast.Applications import NcbiblastnCommandline
import os
import glob
import sys
import re

def get_shared_ids(key_counts, total_db_length):
    '''
    Return the BUSCO ids that are present in all species
    '''
    id_pass = []
    for i in key_counts:
        if key_counts[i] == total_db_length:
            id_pass.append(i)
    return(id_pass)
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


def make_blast_dictionary(results, r, q):
    b_map = {}
    for line in results:
        l = line.strip().split()
        if float(l[2]) >= 30:
            if l[r] not in b_map.keys():
                b_map[l[r]] = [l[q]]
            else:
                b_map[l[r]].append(l[q])
        else:
            pass
    return(b_map)


def run_blast(blast_type, db_path, query_path, outf):
    if blast_type == 'prot':
        if not os.path.exists(db_path + ".nhr"):
            blastdb_cmd = 'makeblastdb -in {0} -dbtype prot'.format(db_path)
            os.system(blastdb_cmd)
        blastp_cline = NcbiblastpCommandline(query=query_path, db=db_path, evalue=0.001, outfmt=6, out=outf)
        stdout, stderr = blastp_cline()
    else:
        if not os.path.exists(db_path + ".nhr"):
            blastdb_cmd = 'makeblastdb -in {0} -dbtype nucl'.format(db_path)
            os.system(blastdb_cmd)
        blastp_cline = NcbiblastnCommandline(query=query_path, db=db_path, evalue=0.001, outfmt=6, out=outf)
        stdout, stderr = blastp_cline()
        
def main(argv):
    if len(argv) != 4:
        sys.stderr.write("Usage: %s <in.sequence.dir> <in.reference.fasta.name> <output.dir> <blast.type>\n" % os.path.basename(sys.argv[0]))
        sys.exit(1)

    db        = sys.argv[1] # directory housing the proteins to get orthologues for
    reference = sys.argv[2] ## in this case whichever one has the fewest proteins
    output_directory = sys.argv[3]
    reference_path = db + "/" + reference ## in this case whichever one has the fewest proteins

    st = db + "/"
    filename, ext = os.path.splitext(reference)
    protein_files = [re.sub(st, '', x) for x in glob.glob(db + "/*" + ext)]
    ## take out the reference sequence
    protein_files.remove(reference)

    all_matches = {}
    key_counts = {}

    for i in protein_files:
        query = i
        query_path = db + "/" + query
        ## make blast dbs if they don't exist already
        run_blast(blast_type = sys.argv[4], db_path = reference_path, query_path = query_path, outf = 'blast_tmp1.tsv')
        run_blast(blast_type = sys.argv[4], db_path = query_path, query_path = reference_path, outf = 'blast_tmp2.tsv')
        result_handle = open("blast_tmp1.tsv")
        b1_map = make_blast_dictionary(result_handle, 1, 0)
        os.remove('blast_tmp1.tsv')
        result_handle = open("blast_tmp2.tsv")
        b2_map = make_blast_dictionary(result_handle, 0, 1)
        os.remove("blast_tmp2.tsv")
        ## return all of the matches
        all_matches[query] = {}
        
        for x in b1_map.keys():
            if x in b2_map.keys():
                if len(b1_map[x]) == 1:
                    if len(b2_map[x]) == 1:
                        if b1_map[x] == b2_map[x]: ## its a reciprocal match!
                            all_matches[query][x] = b2_map[x][0]
                        else:
                            pass ## not reciprocal match
                    else: ## find and rank the reciprocal hit respective to the reference if there are multiple good hits
                    ## only keep things with a rank <= 2
                        for match in b2_map[x]:
                            if x in b2_map:
                                rank = b2_map[x].index(match)
                                if rank <= 2:
                                    all_matches[query][x] = b2_map[x][rank]
                                else:
                                    pass
                            else:
                                pass
                else:
                    match_rank = {}
                    for match in b1_map[x]:
                        ## find the best ranking hit if more than one
                        if match in b2_map[x]:
                            ## create a combined rank to score each hit
                            rank = b2_map[x].index(match) + b1_map[x].index(match) 
                            match_rank[rank] = match
                        else:
                            pass
                    ## keep the hit with the best combined rank
                    vals = list(match_rank.keys())
                    vals.sort()
                    keep_id = match_rank[vals[0]]
                    all_matches[query][x] = keep_id
            else:
                pass                        


    ## for every reference protein sequence count the number of times there are reciprocal hits for each sequence
    key_counts = {}
    for i in all_matches:
        for x in all_matches[i]:
            if x not in key_counts.keys():
                key_counts[x] = 1
            else:
                key_counts[x] += 1
    keeps = get_shared_ids(key_counts = key_counts, total_db_length = len(protein_files))
    ## write out everything
    write_treedb_diction2fasta(matches_dictionary = all_matches,
                                                  fasta_dir = db,
                                                  output_directory = output_directory,
                                                  keep_sequences = keeps)


if __name__ == '__main__':
    main(sys.argv[1:])
