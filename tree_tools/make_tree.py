#!/usr/bin/env python3

import os
import sys
import subprocess
import re
import pysam
from Bio import AlignIO
from Bio import SeqIO

prottest_home = "/home/gotting/src/prottest-3.4.2"

raxml_prot = ['DAYHOFF',
              'DCMUT',
              'JTT',
              'MTREV',
              'WAG',
              'RTREV',
              'CPREV',
              'VT',
              'BLOSUM62',
              'MTMAM',
              'LG',
              'MTART',
              'MTZOA',
              'PMB',
              'HIVB',
              'HIVW',
              'JTTDCMUT',
              'FLU',
              'STMTREV',
              'DUMMY',
              'DUMMY2',
              'AUTO',
              'LG4M',
              'LG4X',
              'PROT_FILE',
              'GTR_UNLINKED',
              'GTR']

def run_prottest(afa, threads, recompute):
    if not os.path.exists('prottest3'):
        os.mkdir('prottest3')
    wd = os.getcwd()
    afa = os.path.join(wd, afa)
    output_file = os.path.join(wd, 'prottest3', os.path.basename(afa))
    if not recompute and os.path.exists(output_file):
        pass
    else:
        command = "./prottest3 -i {} -all-distributions -F -AIC -BIC -AICC -tc 0.5 -o {} -threads {}".format(afa, output_file, threads)
        print('Running:', command)
        p = subprocess.Popen(command, cwd=prottest_home, shell = True)
        os.waitpid(p.pid, 0)
    return(output_file)

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

def run_raxml(model, phylip, threads):
    ''' Runs RaXML using the model, phylip file, and number of threads provided.
    predefined RaXML options used: 
    command1: 
        -f a : rapid Bootstrap analysis and search for best­scoring ML tree in one program run
        -x : Specify an integer number (random seed) and turn on rapid bootstrapping 
             CAUTION:   unlike   in   previous   versions   of   RAxML   will   conduct   rapid   BS  
             replicates under the model of rate heterogeneity you specified via ­m and   
             not by default under CAT 

    command2:
        -f b : draw bipartition information on a tree provided with ­t (typically the bestknown ML tree) 
               based on multiple trees (e.g., from a bootstrap) in a file specified by ­z
    '''
    if not os.path.exists('raxml_trees'):
        os.mkdir('raxml_trees')
    output_dir = os.path.join(os.getcwd(), 'raxml_trees')
    name = re.sub('.phylip', '', os.path.basename(phylip))
    command1 = "raxmlHPC -w {output} -m PROTGAMMAI{model} -n {name} -s {phylip} -f a -x 897543 -p 345232 -N autMRE -T {threads} > {output}/{name}_stdout_raxml1.txt".format(model = model, name = name, phylip = phylip, threads = threads, output = output_dir)
    command2 = "raxmlHPC -w {output} -m PROTGAMMAI{model} -f b -z {output}/RAxML_bootstrap.{name} -t {output}/RAxML_bestTree.{name} -n {name}.BS_tree -T {threads} > {output}/{name}_stdout_raxml2.txt".format(model = model, name = name, threads = threads, output = output_dir)
    print('Running:', command1)
    p = subprocess.Popen(command1, shell = True)
    os.waitpid(p.pid, 0)
    print('Running:', command2)
    p = subprocess.Popen(command2, shell = True)
    os.waitpid(p.pid, 0)

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

    

def main(argv):
    if len(argv) < 2:
        sys.stderr.write("Usage: %s <fastafile> <threads> <recompute>\n" % os.path.basename(sys.argv[0]))
        sys.exit(1)
    fasta   = sys.argv[1]
    threads = sys.argv[2]
    if len(argv) > 2:
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
        rename_fasta_seqs(fasta, out_fasta)
        fasta = out_fasta

    ## align the fasta sequences
    mafft_outfile_name = run_mafft(fasta = fasta, threads = threads, recompute = recompute)

    ## run and parse prottest3 resutls
    prottest_outfile_name = run_prottest(afa = mafft_outfile_name, threads = threads, recompute = recompute)
    model = parse_prottest(prottest_outfile_name)
    ## print out some prottest3 diagnostics
    if model not in raxml_prot:
        exlude_file = open('exluded_alignment_models.txt', 'a')
        print('Model', model, 'not available for raxml, exiting')
        exlude_file.write("{} {}\n".format(mafft_outfile_name, model))
        exlude_file.close()
        sys.exit() # only continue if the model is available in RaXML
    print('Using {} model for RaXML'.format(model))
    include_file = open('included_alignment_models.txt', 'a')
    include_file.write("{} {}\n".format(mafft_outfile_name, model))

    ## convert the alignment fasta file to a phylip
    align           = AlignIO.read(mafft_outfile_name, "fasta")
    phylip_filename = re.sub("afa", "phylip", mafft_outfile_name)
    SeqIO.write(align, phylip_filename, "phylip")

    ## run raxml on the phylip file
    run_raxml(model = model, phylip = phylip_filename, threads = threads)


if __name__ == '__main__':
    main(sys.argv[1:])

