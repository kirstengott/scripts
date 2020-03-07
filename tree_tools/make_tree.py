#!/usr/bin/env python3
import os
import sys
import subprocess
import re
import pysam
from Bio import AlignIO
from Bio import SeqIO
from Bio.Phylo.PAML import yn00
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC, Gapped
import shutil
from itertools import combinations


prottest_home   = "/home/gotting/src/prottest-3.4.2"
jmodeltest_home = "/home/gotting/src/jmodeltest-2.1.10"

## define a dictionary that converts prottest output into raxml matrix models IDs
prottest_raxml = {"JTT" : 'JTT',
                  'LG' : 'LG',
                  'DCMut' : 'DCMUT',
                  'MtREV' : 'MTREV',
                  'MtMam' : 'MTMAM',
                  'MtArt' : 'MTART',
                  'Dayhoff' : 'DAYHOFF',
                  'WAG' : 'WAG',
                  'RtREV' : 'RTREV',
                  'CpREV' : 'CPREV',
                  'Blosum62' : 'BLOSUM62',
                  'VT' : 'VT',
                  'HIVb' : 'HIVB',
                  'HIVw' : 'HIVW',
                  'FLU' : 'FLU'
                  }

stop_codons = ['TAG', 'TAA', 'TGA',
               'tag', 'taa', 'tga']


def check_fasta_headers(fasta):
    '''Check the fasta headers to see if they are longer
       10 characters, if so break and return false'''
    fh = pysam.FastxFile(fasta)
    lengths = set()
    for entry in fh:
        lengths.add(len(entry.name))        
    if max(lengths) > 10:
        return(True)
    else:
        return(False)

def rename_fasta_seqs(fasta, fasta_out):
    '''Renames of fasta file with 'seqN' where N is a sequential number
       and writes out a file with extenstion '_id_map.txt' with mappings
       to the original names fasta sequence names.
       Return: the filename of the id_map for later handling
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
    return(id_map)

def run_prottest(afa, threads, recompute):
    if not os.path.exists('prottest3'):
        os.mkdir('prottest3')
    wd = os.getcwd()
    afa = os.path.join(wd, afa)
    file_base = os.path.basename(afa)
    output_file = os.path.join(wd, 'prottest3', file_base)
    stdout =  os.path.join(wd, 'prottest3', file_base) + "_stdout"
    if not recompute and os.path.exists(output_file):
        pass
    else:
        command = "./prottest3 -i {} -I -F -AIC -BIC -AICC -tc 0.5 -o {} -threads {} > {}".format(afa, output_file, threads, stdout)
        print('Running:', command)
        p = subprocess.Popen(command, cwd=prottest_home, shell = True)
        os.waitpid(p.pid, 0)
    return(output_file)

def parse_prottest(prot_out):
    f = open(prot_out, 'r')
    models = {}
    unique_models = set()
    get_next = 0
    for line in f:
        if 'Best model according to' in line:
            liner = re.sub('Best model according to ', '', line.rstrip()).split()
            score = re.sub(":", "", liner[0])
            model = liner[1]
            unique_models.add(model)
            get_next = 1
        elif get_next == 1:
            ci = line.rstrip().split()[2]
            models[score] = [model, ci]
            get_next = 0
    if len(unique_models) == 1:
        return(list(unique_models)[0])
    else:
        return(models['AICc'][0]) ## if the models don't converge to one, pick the one with best AICc
    
def run_mafft(fasta, threads, recompute):
    if not os.path.exists('mafft'):
        os.mkdir('mafft')
    output = os.path.basename(os.path.splitext(fasta)[0]) + '.afa'
    output_file = os.path.join('mafft', output)
    if not recompute and os.path.exists(output_file):
        pass
    else:
        command = ' '.join(['mafft', '--thread', threads, '--auto', fasta, ">", output_file])
        print('Running:', command)
        p = subprocess.Popen(command, shell=True)
        pid = os.waitpid(p.pid, 0)
        pid = list(pid)[1]
        if pid != 0:
            print('mafft failed, check input fasta file')
            sys.exit()
    return(output_file)

def check_nucleotides_dialign(fasta):
    flag = 0
    fa = pysam.FastxFile(fasta)
    lines = []
    for x in fa:
        seq = ''
        for nuc in x.sequence:
            if nuc not in ['A', 'C', 'T', 'G', 'U']:
                flag = 1
                continue ## remove the bad nucleotides to be treated as gaps by dialign
            else:
                seq += nuc
                
        lines.append(">{}\n{}\n".format(x.name, seq))
    if flag == 1:
        fasta_out = fasta + '.cleaned'
        fh_o = open(fasta_out, 'w')
        [fh_o.write(x) for x in lines]
        return(fasta_out)
    else:
        return(fasta)

    

def run_dialign(fasta, recompute):
    if not os.path.exists('dialigntx'):
        os.mkdir('dialigntx')
    output = os.path.basename(os.path.splitext(fasta)[0]) + '.afa'
    output_file = os.path.join('dialigntx', output)
    if not recompute and os.path.exists(output_file):
        pass
    else:
        dialignconf = '/home/gotting/src/DIALIGN-TX_1.0.2/conf'
        if not os.path.exists(dialignconf):
            print(dialignconf, 'does not exist, fix path')
            sys.exit(1)
        command = "dialign-tx -T {} {} {}".format(dialignconf, fasta, output_file)
        print('Running:', command)
        p = subprocess.Popen(command, shell=True)
        pid = os.waitpid(p.pid, 0)
        pid = list(pid)[1]
        if pid != 0:
            print('dialign failed, check input fasta file')
            sys.exit()
    return(output_file)



def run_trimal(afa, recompute):
    ''' Runs TrimAl to trim the alignments. 
    afa: the output from aligning with mafft
    recompute: whether or not to recopute the trimming
    '''
    output_file = afa + '.trimmed'
    if not recompute and os.path.exists(output_file):
        pass
    else:
        command = "trimal -in {} -out {} -automated1".format(afa, output_file)
        print('Running:', command)
        p = subprocess.Popen(command, shell=True)
        pid = os.waitpid(p.pid, 0)
        pid = list(pid)[1]
        if pid != 0:
            print('trimal failed, check input fasta file')
            sys.exit()
    return(output_file)
        

def run_yn00(phylip, recompute, id_map = False):
    if not os.path.exists('yn_out'):
        os.mkdir('yn_out')
    path_base = re.sub('.phylip', '', os.path.basename(phylip))
    output_filename = os.path.join(os.getcwd(), 'yn_out', path_base)
    infile = os.path.join(os.getcwd(), phylip)
    if not recompute and os.path.exists(output_filename):
        pass
    else:
        #ctl_file = output_filename + '_CTL'
        yn = yn00.Yn00()
        yn.alignment = infile
        yn.out_file = output_filename + '_full'
        yn.working_dir = os.getcwd()
        try:
            yn.run(verbose = True, command = '/usr/local/bin/yn00')
            results = yn00.read(output_filename + "_full")
            #print(results.keys())
        except:
            print('Error parsing results for', infile, 'check', output_filename)
            results = False
        if results:
            if id_map:
                ids = {}
                with open(id_map, 'r') as f:
                    for line in f:
                        (key, val) = line.strip().split()
                        ids[key] = val
            output_file = open(output_filename, 'w')
            output_file.write("seq1,seq2,dn,dn_se,ds,ds_se\n")
            for i in combinations(results.keys(), 2):
                seq1 = i[0]
                seq2 = i[1]
                res = results[seq1][seq2]['YN00']
                if id_map:
                    outputs = [ids[seq1], ids[seq2], res['dN'], res['dN SE'], res['dS'], res['dS SE']]

                else:
                    outputs = [seq1, seq2, res['dN'], res['dN SE'], res['dS'], res['dS SE']]
                outputs = ",".join([str(x) for x in outputs])
                output_file.write(outputs + "\n")



def run_raxml(model, phylip, threads, recompute, cds = False):
    ''' Runs RaXML using the model, phylip file, and number of threads provided.
    predefined RaXML options used: 
    command1: 
        -f a : rapid Bootstrap analysis and search for best­scoring ML tree in one program run
        -x : Specify an integer number (random seed) and turn on rapid bootstrapping 
             CAUTION:   unlike   in   previous   versions   of   RAxML   will   conduct   rapid   BS  
             replicates under the model of rate heterogeneity you specified via ­m and   
             not by default under CAT 
    ## not used
    command2:
        -f b : draw bipartition information on a tree provided with ­t (typically the bestknown ML tree) 
               based on multiple trees (e.g., from a bootstrap) in a file specified by ­z
    '''
    if not os.path.exists('raxml_trees'):
        os.mkdir('raxml_trees')
    output_dir = os.path.join(os.getcwd(), 'raxml_trees')
    name = re.sub('.phylip', '', os.path.basename(phylip))
    output_filename = 'RAxML_bipartitionsBranchLabels.' + name
    output_filename = os.path.join(output_dir, output_filename)
    if not recompute and os.path.exists(output_filename):
        pass
    else:
        if cds:
             command1 = "raxmlHPC -w {output} -m {model} -n {name} -s {phylip} -f a -x 897543 -p 345232 -N autoMRE -T {threads} > {output}/{name}_stdout_raxml1.txt".format(model = model, name = name, phylip = phylip, threads = threads, output = output_dir)
        else:
            command1 = "raxmlHPC -w {output} -m PROTGAMMA{model} -n {name} -s {phylip} -f a -x 897543 -p 345232 -N autoMRE -T {threads} > {output}/{name}_stdout_raxml1.txt".format(model = model, name = name, phylip = phylip, threads = threads, output = output_dir)
        print('raxml command:', command1)
        p = subprocess.Popen(command1, shell = True)
        os.waitpid(p.pid, 0)
    return(output_filename)


def run_fasttree(phylip, threads, recompute, cds = False):
    ''' Runs FastTree, phylip file, and number of threads provid.
    '''
    if not os.path.exists('fastree'):
        os.mkdir('fastree')
    output_dir = os.path.join(os.getcwd(), 'fastree')
    name = re.sub('.phylip', '', os.path.basename(phylip))
    output_filename = 'Fastree.' + name
    output_filename = os.path.join(output_dir, output_filename)
    if not recompute and os.path.exists(output_filename):
        pass
    else:
        if cds == True:
            command1 = "fasttree -gtr -nt < {phylip} > {output} ".format(phylip = phylip, output = output_filename)
        else:
            command1 = "fasttree < {phylip} > {output} ".format(phylip = phylip, output = output_filename)
        print('fastree command:', command1)
        p = subprocess.Popen(command1, shell = True)
        os.waitpid(p.pid, 0)
    return(output_filename)


def rename_newick(id_map, nwk_in):
    print(nwk_in)
    nwk_f = open(nwk_in, 'r')
    id_map_f = open(id_map, 'r')
    nwk = nwk_f.readlines()[0]
    nwk_f.close()
    for line in id_map_f:
        line = line.rstrip().split()
        temp_id = line[0] + ':'
        new_id = re.sub(",", ".", line[1]) + ":"
        nwk = re.sub(temp_id, new_id, nwk)
    id_map_f.close()
    if not os.path.exists('final_trees'):
        os.mkdir('final_trees')
    nwk_out = os.path.join('final_trees', os.path.basename(nwk_in))
    nwk_out_h = open(nwk_out, 'w')
    nwk_out_h.write(nwk)
    nwk_out_h.close()
    return(nwk_out)

def clean_phylip_yn00(align_outfile_name, recompute):
    stop_flag = False
    phylip_yn00_filename = re.sub("afa", "phylip.yn00", align_outfile_name)
    if recompute == True and os.path.exists(phylip_yn00_filename):
        os.remove(phylip_yn00_filename)
    if not os.path.exists(phylip_yn00_filename):
        align           = AlignIO.read(align_outfile_name, "fasta")
        for record in align:
            len_r = len(record.seq)
            last_codon = record.seq[(len_r-3):len_r]
            ## remove the stop codon for yn00
            if last_codon in stop_codons:
                stop_flag = True
        if stop_flag == True:
            for record in align:
                record.seq = record.seq[:-3]
        SeqIO.write(align, phylip_yn00_filename, "phylip-sequential")
    return(phylip_yn00_filename)
        

def clean_to_phylip(align_outfile_name, recompute):
        ## convert to phylip and remove redundant  sequences
    phylip_filename = re.sub("afa.trimmed", "phylip", align_outfile_name)
    if recompute == True and os.path.exists(phylip_filename):
        os.remove(phylip_filename)
    if not os.path.exists(phylip_filename):
        try:
            align           = AlignIO.read(align_outfile_name, "fasta")
            uniq_seqs = {} ## sequence is the key, values are the names
            ## make the alignment only have unique sequences
            for record in align:
                seq = str(record.seq)
                sid = str(record.id)
                ## remove empty records that didn't align
                if set(seq) == set("-"):
                    continue
                if not seq in uniq_seqs:
                    uniq_seqs[seq] = sid
                else:
                    add_record = "," + sid
                    uniq_seqs[seq] += add_record
            align1 = AlignIO.MultipleSeqAlignment([], Gapped(IUPAC.unambiguous_dna, "-"))
            for x in uniq_seqs:
                align1.add_sequence(uniq_seqs[x], x)
            if align1.get_alignment_length() >= 30:
                SeqIO.write(align, phylip_filename, "phylip-sequential")
        except:
            print(align_outfile_name, 'too short to continue after trimming')
    return(phylip_filename)

def main_cds_tree(fasta, threads, recompute, id_map_filename = False):
    fasta                = check_nucleotides_dialign(fasta)
    dialign_outfile_name = run_dialign(fasta = fasta, recompute = recompute)
    align_outfile_name   = run_trimal(afa = dialign_outfile_name, recompute = recompute)
    phylip_filename      = clean_to_phylip(align_outfile_name = align_outfile_name, recompute = recompute)
            
    if '-yn' in sys.argv:
        phylip_yn00_filename = clean_phylip_yn00(align_outfile_name = dialign_outfile_name, recompute = recompute)
        run_yn00(phylip = phylip_yn00_filename, recompute = recompute, id_map = id_map_filename)
    if '-a' in sys.argv:
        sys.exit()
        
    ## align the phylip file
    if '-fasttree' in sys.argv:
        tree_filename = run_fasttree(phylip = phylip_filename, threads = threads, recompute = recompute, cds = True)
    else:
        model = 'GTRGAMMA'
        tree_filename = run_raxml(model = model, phylip = phylip_filename, threads = threads, recompute = recompute, cds = True)
    return(tree_filename)

    
def main_protein_tree(fasta, threads, recompute):
    mafft_outfile_name = run_mafft(fasta = fasta, threads = threads, recompute = recompute)
    align_outfile_name = run_trimal(afa = mafft_outfile_name, recompute = recompute)
    phylip_filename    = clean_to_phylip(align_outfile_name = align_outfile_name, recompute = recompute)
    print(phylip_filename)
    if '-a' in sys.argv:
        sys.exit()
        
    ## align the phylip file
    if '-fasttree' in sys.argv:
        tree_filename = run_fasttree(phylip = phylip_filename, threads = threads, recompute = recompute)
    else:
        ## run and parse prottest3 resutls
        prottest_outfile_name = run_prottest(afa = phylip_filename, threads = threads, recompute = recompute)
        if '-p' in sys.argv:
            sys.exit()
        model_init = parse_prottest(prottest_outfile_name)
        model_all = model_init.split("+")

        ## check to see if the model is available to raxml, otherwise report it
        try:
            if len(model_all) > 1:
                if model_all[1] == 'I':
                    model = 'I' + prottest_raxml[model_all[0]]
                else:
                    model = prottest_raxml[model_all[0]]
                    print('unkown model:', model_all, 'using', model[0])
            else:
                model = prottest_raxml[model_all[0]]
            print('Using {} model for RaXML'.format(model))
            include_file = open('included_alignment_models.txt', 'a')
            include_file.write("{} {}\n".format(mafft_outfile_name, model))
        except:
            exlude_file = open('exluded_alignment_models.txt', 'a')
            print('Model', model_init, 'not available for raxml, having raxml estimate')
            exlude_file.write("{} {}\n".format(phylip_filename, model_init))
            exlude_file.close()
            sys.exit()
        tree_filename = run_raxml(model = model, phylip = phylip_filename, threads = threads, recompute = recompute)
    return(tree_filename)
        


def main(argv):
    if len(argv) < 2:
        sys.stderr.write('''Usage: %s <fastafile> <threads> <-r>
        <fastafile> : the fasta file of orthologs to make a tree of
        <threads> : the number of threads to use
        <-r> : recompute everything
        <-a> : stop after aligning
        <-p> : stop after prottest
        <-fasttree>: use fasttree instead of raxml/prottest
        <-cds>: input type is cds (protein by default)
        <-yn> : run yn00 to get dn/ds for cds sequences\n''' % os.path.basename(sys.argv[0]))
        sys.exit(1)
    fasta   = sys.argv[1]
    threads = sys.argv[2]
    if '-r' in sys.argv:
        recompute = True
    else:
        recompute = False
        
    ## first check the fasta file for names lengths
    fa_header_flag = check_fasta_headers(fasta)
    ## if the headers are too long, make new files with shorter headers
    if fa_header_flag:
        renamed_fa_dir = os.path.dirname(fasta) + '_renamed'
        if not os.path.exists(renamed_fa_dir):
            os.mkdir(renamed_fa_dir)
        out_fasta = os.path.join(renamed_fa_dir, os.path.basename(fasta))
        id_map_filename = rename_fasta_seqs(fasta, out_fasta)
        fasta = out_fasta
    else:
        id_map_filename = False

    ## implement different subroutines for cds/protein sequences
    if '-cds' in sys.argv:    
        tree_filename = main_cds_tree(fasta = fasta, threads = threads, recompute = recompute, id_map_filename = id_map_filename)
    else:
        tree_filename = main_protein_tree(fasta = fasta, threads = threads, recompute = recompute)
        
    ## pull the tree in and fix the newick back to names from the original fasta file
    print(tree_filename)
    if id_map_filename:
        nwk_out = rename_newick(id_map = id_map_filename, nwk_in = tree_filename)
        print('Final tree:', nwk_out)
    else:
        nwk_out = os.path.join('final_trees', os.path.basename(tree_filename))
        shutil.copyfile(tree_filename, nwk_out)
        print('Final tree:', nwk_out)

if __name__ == '__main__':
    main(sys.argv[1:])

