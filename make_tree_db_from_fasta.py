#!/usr/bin/env python
from Bio.Blast.Applications import NcbiblastpCommandline
import os
import glob
import make_busco_tree_db
import sys
import re

def make_blast_dictionary(results, r, q):
    b_map = {}
    for line in results:
        l = line.strip().split()
        if float(l[2]) > 70:
            b_map[l[r]] = l[q]
        else:
            pass
    return(b_map)

def main(argv):
    if len(argv) != 3:
        sys.stderr.write("Usage: %s <in.protein.dir> <in.reference.fasta.name> <output.dir>\n" % os.path.basename(sys.argv[0]))
        sys.exit(1)

    db        = sys.argv[1] # directory housing the proteins to get orthologues for
    reference = sys.argv[2] ## in this case whichever one has the fewest proteins
    output_directory = sys.argv[3]



    reference_path = db + "/" + reference ## in this case whichever one has the fewest proteins
    st = db + "/"
    protein_files = [re.sub(st, '', x) for x in glob.glob(db + "/*fa")]                 
    ## take out the reference sequence
    protein_files.remove(reference)

    all_matches = {}
    key_counts = {}

    for i in protein_files:
        query = i
        query_path = db + "/" + query
        ## make blast dbs if they don't exist already
        if not os.path.exists(reference_path + ".phr"):
            blastdb_cmd = 'makeblastdb -in {0} -dbtype prot'.format(reference_path)
            os.system(blastdb_cmd)
        if not os.path.exists(query_path + ".phr"):
            blastdb_cmd = 'makeblastdb -in {0} -dbtype prot'.format(query_path)
            os.system(blastdb_cmd)
        ## run the blasts    
        blastp_cline = NcbiblastpCommandline(query=query_path, db=reference_path, evalue=0.001, outfmt=6, out="blast_tmp1.tsv")
        stdout, stderr = blastp_cline()
        result_handle = open("blast_tmp1.tsv")
        b1_map = make_blast_dictionary(result_handle, 1, 0)
        os.remove('blast_tmp1.tsv')
        blastp_cline = NcbiblastpCommandline(query=reference_path, db=query_path, evalue=0.001, outfmt=6, out="blast_tmp2.tsv")
        stdout, stderr = blastp_cline()
        result_handle = open("blast_tmp2.tsv")
        b2_map = make_blast_dictionary(result_handle, 0, 1)
        os.remove("blast_tmp2.tsv")
        ## return all of the matches
        all_matches[query] = {}
        for x in b1_map.keys():
            if b1_map[x] == b2_map[x]: ## its a reciprocal match!
                all_matches[query][x] = b2_map[x]
                if x in key_counts:
                    key_counts[x] += 1
                else:
                    key_counts[x] = 1
            else:
                pass
            
    keeps = make_busco_tree_db.get_shared_ids(key_counts = key_counts, total_db_length = len(protein_files))

    ## write out everything
    make_busco_tree_db.write_treedb_diction2fasta(matches_dictionary = all_matches,
                                                  fasta_dir = db,
                                                  output_directory = output_directory,
                                                  keep_sequences = keeps)


    

if __name__ == '__main__':
    main(sys.argv[1:])
