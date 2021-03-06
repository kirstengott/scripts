#!/usr/bin/env python3

import os, sys, re, subprocess
from Bio import AlignIO, SeqIO
import pysam
from Bio.Alphabet import IUPAC, Gapped

def parse_modeltest(in_file):
    in_file = in_file
    model = ''
    with open(in_file, 'r') as fh:
        for line in fh:
            line = line.strip()
            if 'Best-fit model:' in line:
                    model = line.split(" ")[2]
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
    #output_filename = name + '.raxml.bestTree'
    #output_filename = os.path.join(output_dir, output_filename)
    output_filename = 'RAxML_bipartitionsBranchLabels.' + name
    command1 = "raxmlHPC-PTHREADS-AVX -m GTRGAMMA -n {outp} -s {phylip} -f a -x 897543 -N autoMRE -p 345232 -T {threads}".format(phylip = phylip, threads = threads, outp = name)
    #command1 = "raxml-ng --all --msa {phylip} --model {model} --prefix {outp} --seed 2 --threads {threads} --brlen scaled > {outp}_stdout_raxml1.txt".format(model = model, phylip = phylip, threads = threads, outp = prefix)
    print('raxml command:', command1)
    p = subprocess.Popen(command1, shell = True)
    os.waitpid(p.pid, 0)
    return(output_filename)

def run_iqtree(phylip, threads, partition, cds = False):
    ''' Runs iqtree with phylip file, and number of threads provided.
    '''
    rax_out = 'combo_iqtree_tree'
    if not os.path.exists(rax_out):
        os.mkdir(rax_out)
    output_dir = os.path.join(os.getcwd(), rax_out)
    name = re.sub('.phylip', '', os.path.basename(phylip))
    prefix = output_dir + "/" + name
    output_filename = name + '.treefile'
    output_filename = os.path.join(output_dir, output_filename)

    command1 = "iqtree2 -s {phylip} -nt {threads} -bb 1000 -alrt 1000 -pre {outp} -p {part_file}> {outp}_stdout_iqtree.txt".format(phylip = phylip, threads = threads, outp = prefix, part_file = partition)
    print('iqtree command:', command1)
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
    nwk_out = os.path.join('combo_iqtree_tree', os.path.basename(nwk_in))
    nwk_out_h = open(nwk_out, 'w')
    nwk_out_h.write(nwk)
    nwk_out_h.close()
    return(nwk_out)

def main(argv):
    if len(argv) < 2:
        sys.stderr.write('''Usage: %s <fastafile> <threads> <-r>
        <aln_directory>  directory housing multiple alignment output
        <threads>        the number of threads to use
        <seq_id_dir>     directory with seq ids if applicable\n''' % os.path.basename(sys.argv[0]))
        sys.exit(1)

    fasta_dir       = sys.argv[1]
    threads         = sys.argv[2]
    if len(argv) == 3:
        id_maps         = sys.argv[3]
    else:
        id_maps = False
    
    fastas    = [x for x in os.listdir(fasta_dir) if 'afa.trimmed' in x]

    seqs_all   = {}
    coords_all = {}
    start = 1
    total_len = 0

    sp_name = set()
    gene_lens = {}

    if id_maps:
        for fasta in fastas:
            gene_name_f = os.path.join(id_maps, re.sub('.afa.trimmed', '', fasta) + "_id_map.txt")
            gene_name_dict = {}
            with open(gene_name_f, 'r') as fh:
                for line in fh:
                    line = line.strip().split()
                    gene_name_dict[line[0]] = line[1]
            fasta = os.path.join(fasta_dir, fasta)
            fa = pysam.FastxFile(fasta)

            num_seq = 0
            for x in fa:
                num_seq += 1
                name = gene_name_dict[x.name]
                sp_name.add(name)
                if fasta not in seqs_all:
                    seqs_all[fasta] = {}

                if name not in seqs_all[fasta]:
                    seqs_all[fasta][name] = x.sequence

            fa.close()
    else:
        for fasta in fastas:
            fasta = os.path.join(fasta_dir, fasta)
            fa = pysam.FastxFile(fasta)
            num_seq = 0
            for x in fa:
                num_seq += 1
                name = x.name
                sp_name.add(name)
                if fasta not in seqs_all:
                    seqs_all[fasta] = {}

                if name not in seqs_all[fasta]:
                    seqs_all[fasta][name] = x.sequence

            fa.close()

        
        
    for i in seqs_all:
        coord_flag = 0
        len_1 = False
        for x in seqs_all[i]: ## pull out the length of the gene
            gene_len = int(len(seqs_all[i][x]))
            if start == 1:
                end = gene_len
            else:
                end = start + gene_len -1
            coords_all[i] = [start, end]
            start = end + 1 
            total_len += gene_len
            break
        for species in sp_name: ## add in missing sequences with gaps
            if not species in seqs_all[i]:
                seq = "-" * gene_len
                seqs_all[i][species] = seq


    seqs_all_t = {}
    for i in seqs_all:
        for x in seqs_all[i]:
            if x not in seqs_all_t:
                seqs_all_t[x] = seqs_all[i][x]
            else:
                seqs_all_t[x] += seqs_all[i][x]

    
    #print(seqs_all_t['ICBG712'])
    #print(total_len)
    #print(len(seqs_all_t['ICBG712']))
    #sys.exit()

    partition_f = 'concatenated_alignments.part'
    partition = open(partition_f, 'w')
    partition.write("#nexus\nbegin sets;\n")
    models = []
    ## make the nexus file for partitioned tree

    for i in coords_all:
        gene_name = re.sub('.afa.trimmed', '', os.path.basename(i))
        modeltest_out = 'iqtrees/' + gene_name + ".log"
        model = parse_modeltest(modeltest_out)
        outstring = "\tcharset {gene} = {start}-{end};\n".format(gene = gene_name,
                                               start = coords_all[i][0],
                                               end = coords_all[i][1])
        partition.write(outstring)
        models.append(model + ":" + gene_name)
    all_models = ",".join(models)
    partition.write("\tcharpartition mine = " + all_models + ';' + "\nend;\n")
    partition.close()
    
    
    align = AlignIO.MultipleSeqAlignment([], Gapped(IUPAC.unambiguous_dna, "-"))
    phylip_filename = 'concatenated_alignments.phylip'
    fasta_filename  = 'concatenated_alignments.fasta'
    count = 0
    id_map = 'concatenated_alignments_id_map.txt'
    ids = open(id_map, 'w')

    

    
    
    for x in seqs_all_t: 
        seqname = "seq" + str(count)
        ids.write(seqname + "\t" + x + "\n")
        align.add_sequence(seqname, seqs_all_t[x])
        count += 1
    ids.close()

    SeqIO.write(align, phylip_filename, "phylip-sequential")
    SeqIO.write(align, fasta_filename, "fasta")

    tree_filename = run_iqtree(phylip = phylip_filename, threads = 'AUTO', partition = partition_f, cds = False)
    ## pull the tree in and fix the newick back to names from the original fasta file
    if not os.path.exists('combo_iqtree_tree'):
        os.mkdir('combo_iqtree_tree')
    if id_map:
        nwk_out = rename_newick(id_map = id_map, nwk_in = tree_filename)
    else:
        nwk_out = os.path.join('combo_raxml_tree', os.path.basename(tree_filename))
        shutil.copyfile(tree_filename, nwk_out)
    print('Final tree:', nwk_out)


if __name__ == '__main__':
    main(sys.argv[1:])
