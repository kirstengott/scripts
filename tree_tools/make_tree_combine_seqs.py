#!/usr/bin/env python3

import os, sys, re, subprocess
from Bio import AlignIO, SeqIO
import pysam
from Bio.Alphabet import IUPAC, Gapped

def parse_modeltest(in_file):
    in_file = in_file + '.out'
    next_line = False
    model = ''
    with open(in_file, 'r') as fh:
        for line in fh:
            line = line.strip()
            if 'Best model according to AICc' in line:
                next_line = True
                continue
            elif next_line == True:
                if '---' in line:
                    continue
                else:
                    model = re.sub('Model:', '', line)
                    model = model.strip()
                    break
            else:
                continue
    return(model)

def run_raxml(model, phylip, threads):
    ''' Runs RaXML using the model, phylip file, and number of threads provided.
    '''
    rax_out = 'combo_raxml_tree'
    if not os.path.exists(rax_out):
        os.mkdir(rax_out)
    output_dir = os.path.join(os.getcwd(), rax_out)
    name = re.sub('.phylip', '', os.path.basename(phylip))
    output_filename = name + '.raxml.bestTree'
    output_filename = os.path.join(output_dir, output_filename)
    prefix = output_dir + "/" + name
    command1 = "raxml-ng --all --msa {phylip} --model {model} --prefix {outp} --seed 2 --threads {threads} --brlen scaled > {outp}_stdout_raxml1.txt".format(model = model, phylip = phylip, threads = threads, outp = prefix)
    print('raxml command:', command1)
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
    nwk_out = os.path.join('combo_raxml_tree', os.path.basename(nwk_in))
    nwk_out_h = open(nwk_out, 'w')
    nwk_out_h.write(nwk)
    nwk_out_h.close()
    return(nwk_out)

def main(argv):
    if len(argv) < 2:
        sys.stderr.write('''Usage: %s <fastafile> <threads> <-r>
        <aln_directory>  directory housing multiple alignment output
        <threads>        the number of threads to use
        <id_map>         id map filename, or 'False'\n''' % os.path.basename(sys.argv[0]))
        sys.exit(1)

    fasta_dir       = sys.argv[1]
    id_map_filename = sys.argv[3]
    threads         = sys.argv[2]
    fastas    = [x for x in os.listdir(fasta_dir) if 'afa.trimmed' in x]

    seqs_all   = {}
    coords_all = {}
    start = 1
    total_len = 0
    for fasta in fastas:
        fasta = os.path.join(fasta_dir, fasta)
        fa = pysam.FastxFile(fasta)
        coord_flag = 0
        for x in fa:
            if x.name in seqs_all:
                seqs_all[x.name] += x.sequence
            else:
                seqs_all[x.name] = x.sequence
            if coord_flag == 0:
                if start == 1:
                    end = len(x.sequence)
                else:
                    end = start + len(x.sequence) -1
                    
                coords_all[fasta] = [start, end]
                start = end + 1 
                total_len += len(x.sequence)
                coord_flag = 1
        fa.close()

   # print(total_len)
   # print(coords_all['dialigntx/133766at4890.fa.afa.trimmed'])
   # print(coords_all['dialigntx/133766at4890.fa.afa.trimmed'][1] - coords_all['dialigntx/133766at4890.fa.afa.trimmed'][0] + 1, 'should be 2376')

    partition_f = 'concatenated_alignments.part'
    partition = open(partition_f, 'w')
    for i in coords_all:
        gene_name = re.sub('.afa.trimmed', '', os.path.basename(i))
        modeltest_out = 'modeltest/' + os.path.basename(i) + ".COMPLETE"
        model = parse_modeltest(modeltest_out)
        outstring = "{model}, {gene}={start}-{end}\n".format(model = model, gene = gene_name,
                                               start = coords_all[i][0],
                                               end = coords_all[i][1])
        partition.write(outstring)
    partition.close()
    
    align = AlignIO.MultipleSeqAlignment([], Gapped(IUPAC.unambiguous_dna, "-"))
    phylip_filename = 'concatenated_alignments.phylip'
    fasta_filename  = 'concatenated_alignments.fasta'
    for x in seqs_all:
        align.add_sequence(x, seqs_all[x])
    SeqIO.write(align, phylip_filename, "phylip-sequential")
    SeqIO.write(align, fasta_filename, "fasta")
    tree_filename = run_raxml(model = partition_f,
                               phylip = phylip_filename,
                               threads = threads)

    ## pull the tree in and fix the newick back to names from the original fasta file
    if not os.path.exists('combo_raxml_tree'):
        os.mkdir('combo_raxml_tree')
    if id_map_filename:
        nwk_out = rename_newick(id_map = id_map_filename, nwk_in = tree_filename)
    else:
        nwk_out = os.path.join('combo_raxml_tree', os.path.basename(tree_filename))
        shutil.copyfile(tree_filename, nwk_out)
    print('Final tree:', nwk_out)


if __name__ == '__main__':
    main(sys.argv[1:])
