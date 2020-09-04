#!/usr/bin/env python3
import os, sys, subprocess, re
import pysam
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO
from Bio import SeqIO
from Bio.Phylo.PAML import yn00
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC, Gapped, generic_dna
import shutil
from itertools import combinations
from pathlib import Path

script_dir = os.path.dirname(os.path.realpath(__file__))


stop_codons = ['TAG', 'TAA', 'TGA',
               'tag', 'taa', 'tga']

dialignconf = os.path.join(script_dir, 'conf/DIALIGN-TX_1.0.2_conf')


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


def run_modeltest(in_file, output_file, threads, seq_type, recompute = False):
    output_here = output_file + '.out'
    if not recompute and os.path.exists(output_here):
        pass
    else:
        if not os.path.exists('modeltest'):
            os.mkdir('modeltest')
        command = "modeltest-ng -i {input_f} -d {stype} -o {out} -p {threads} -T raxml".format(input_f = in_file, threads = threads, out = output_file, stype = seq_type)
        print('Running:', command, "to make", output_file)
        os.system(command)

def parse_modeltest(in_file):
    in_file = in_file + '.out'
    next_line = False
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
            
    
def run_phyml(phylip, model, recompute):
    outfile = phylip + '_phyml_tree.txt'
    return_file = os.path.join('phyml_trees', os.path.basename(outfile))
    if not os.path.exists('phyml_trees'):
        os.mkdir('phyml_trees')
    if not recompute and os.path.exists(return_file):
        pass
    else:
        command = "phyml -i {phylip} -q -d nt -m {model} -a e -b 100".format(phylip = phylip, model = model)
        print('Running:', command)
        p = subprocess.Popen(command, shell = True)
        os.waitpid(p.pid, 0)
        shutil.copyfile(outfile, return_file)
    return(return_file)
    


    
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
        command = "dialign-tx -T {} {} {}".format(dialignconf, fasta, output_file)
        print('Running:', command)
        p = subprocess.Popen(command, shell=True)
        pid = os.waitpid(p.pid, 0)
        pid = list(pid)[1]
        if pid != 0:
            print('dialign failed, check input fasta file')
            sys.exit()
    return(output_file)


def run_prank(fasta, recompute):
    if not os.path.exists('prank'):
        os.mkdir('prank')
    output_prefix = os.path.join('prank', os.path.splitext(os.path.basename(os.path.splitext(fasta)[0]))[0])
    output_file   = output_prefix + '.best.fas'
    if not recompute and os.path.exists(output_file):
        pass
    else:
        command = "prank -d={} -o={} -showanc -showevents -codon -F".format(fasta, output_prefix)
        print('Running:', command)
        p = subprocess.Popen(command, shell=True)
        pid = os.waitpid(p.pid, 0)
        pid = list(pid)[1]
        if pid != 0:
            print('prank failed, check input fasta file')
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
        yn.run(verbose = True)
        results = yn00.read(output_filename + "_full")
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
        else:
            print('Error parsing results for', infile, 'check', output_filename)
            results = False



def run_raxml(model, phylip, threads, recompute, cds = False):
    ''' Runs RaXML using the model, phylip file, and number of threads provided.
    '''
    if not os.path.exists('raxml_trees'):
        os.mkdir('raxml_trees')
    output_dir = os.path.join(os.getcwd(), 'raxml_trees')
    name = re.sub('.phylip', '', os.path.basename(phylip))
    prefix = output_dir + "/" + name
    output_filename = prefix + ".raxml.bestTree"
    if not recompute and os.path.exists(prefix + '.raxml.bootstraps'):
        pass
    else:
        command1 = "raxml-ng --bootstrap --msa {phylip} --model {model} --prefix {outp} --seed 2 --bs-trees 1000 --threads {threads} > {outp}_stdout_raxml1.txt".format(model = model, phylip = phylip, threads = threads, outp = prefix)
        print('raxml command:', command1)
        p = subprocess.Popen(command1, shell = True)
        os.waitpid(p.pid, 0)
    return(output_filename)






def run_iqtree(phylip, threads, recompute, cds = False):
    ''' Runs iqtree with phylip file, and number of threads provided.
    '''
    if not os.path.exists('iqtrees'):
        os.mkdir('iqtrees')
    output_dir = os.path.join(os.getcwd(), 'iqtrees')
    name = re.sub('.phylip', '', os.path.basename(phylip))
    prefix = output_dir + "/" + name
    output_filename = name + '.treefile'
    output_filename = os.path.join(output_dir, output_filename)
    if not recompute and os.path.exists(output_filename):
        pass
    else:
        command1 = "iqtree -s {phylip} -nt {threads} -bb 1000 -alrt 1000 -pre {outp} > {outp}_stdout_iqtree.txt".format(phylip = phylip, threads = threads, outp = prefix)
        print('iqtree command:', command1)
        p = subprocess.Popen(command1, shell = True)
        os.waitpid(p.pid, 0)
    return(output_filename)






def run_fasttree(phylip, threads, recompute, cds = False):
    ''' Runs FastTree, phylip file, and number of threads provid.
    '''
    if not os.path.exists('fastree'):
        os.mkdir('fastree')
    output_dir = os.path.join(os.getcwd(), 'fastree')
    name = re.sub("\.*$", '', (os.path.basename(fasta)))
    output_filename = name
    output_filename = os.path.join(output_dir, output_filename)
    if not recompute and os.path.exists(output_filename):
        pass
    else:
        if cds == True:
            command1 = "fasttree -gtr -nt < {phylip} > {output} ".format(phylip = fasta, output = output_filename)
        else:
            command1 = "fasttree < {phylip} > {output} ".format(phylip = fasta, output = output_filename)
        print('fastree command:', command1)
        p = subprocess.Popen(command1, shell = True)
        os.waitpid(p.pid, 0)
    return(output_filename)


def rename_newick(id_map, nwk_in):
    #print(nwk_in)
    nwk_f = open(nwk_in, 'r')
    id_map_f = open(id_map, 'r')
    nwk = nwk_f.readlines()[0]
    nwk_f.close()
    for line in id_map_f:
        line = line.rstrip().split()
        temp_id = line[0] + ':'
        if len(re.findall(temp_id, line[1])):
                print('Mulitple taxa placement:', new_id)
                continue
        else:
            new_id = re.sub(",", ".", line[1]) + ":"
            nwk = re.sub(temp_id, new_id, nwk)
    id_map_f.close()
    nwk_out = os.path.join('final_trees', os.path.basename(nwk_in))
    nwk_out_h = open(nwk_out, 'w')
    nwk_out_h.write(nwk)
    nwk_out_h.close()
    return(nwk_out)




def clean_fasta_codons(fasta, recompute):
    stop_flag = False
    fasta_clean = fasta + '.codons'
    if recompute == True and os.path.exists(fasta_clean):
        os.remove(fasta_clean)
    if not os.path.exists(fasta_clean):
        fa_in = pysam.FastxFile(fasta)
        ## find stop codons in the last position
        stop_indices = []
        codon_3_flag = False
        for record in fa_in:
            if(len(record.sequence)%3 != 0):
                codon_3_flag = True
            codons = re.findall(r"(.{3})", str(record.sequence))
            index = 0
            for x in codons:
                if x in stop_codons:
                    if index not in stop_indices:
                        stop_indices.append(index)
                    index += 1
                else:
                    index += 1
                    continue
        fa_in = pysam.FastxFile(fasta)
        if len(stop_indices) > 0 or codon_3_flag == True:
            fa_out = open(fasta_clean, 'w')
            for record in fa_in:
                codons = re.findall(r"(.{3})", str(record.sequence))
                index = 0
                new_seq = ''
                for x in codons:
                    if index in stop_indices:
                        index += 1
                        continue
                    else:
                        index += 1
                        new_seq += x
                temp_seq = ">{}\n{}\n".format(record.name, new_seq)
                fa_out.write(temp_seq)
            fa_out.close()
        else:
            fasta_clean = fasta

    return(fasta_clean)






def clean_phylip_yn00(align_outfile_name, recompute):
    stop_flag = False
    phylip_yn00_filename = re.sub("afa", "phylip.yn00", align_outfile_name)
    if recompute == True and os.path.exists(phylip_yn00_filename):
        os.remove(phylip_yn00_filename)
    if not os.path.exists(phylip_yn00_filename):
        align           = AlignIO.read(align_outfile_name, "fasta")
        ## find stop codons in the last position
        stop_indices = []
        codon_3_flag = False
        for record in align:
            if(len(record.seq)%3 != 0):
                codon_3_flag = True
            codons = re.findall(r"(.{3})", str(record.seq))
            index = 0
            for x in codons:
                if x in stop_codons:
                    if index not in stop_indices:
                        stop_indices.append(index)
                    index += 1
                else:
                    index += 1
                    continue
        if len(stop_indices) > 0 or codon_3_flag == True:
            align1 = AlignIO.MultipleSeqAlignment([], Gapped(IUPAC.unambiguous_dna, "-"))
            for record in align:
                codons = re.findall(r"(.{3})", str(record.seq))
                index = 0
                new_seq = ''
                for x in codons:
                    if index in stop_indices:
                        index += 1
                        continue
                    else:
                        index += 1
                        new_seq += x
                temp_seq = SeqRecord(Seq(new_seq, generic_dna), id = record.id)
                align1.append(temp_seq)
            AlignIO.write(align1, phylip_yn00_filename, "phylip-sequential")
        else:
            AlignIO.write(align, phylip_yn00_filename, "phylip-sequential")
    return(phylip_yn00_filename)
        



def clean_dialign_fa(fasta, recompute):
        ## convert to phylip and remove redundant  sequences
    lower_nucs = ['a', 't', 'g', 'c']
    out_filename = fasta + ".clean"
    if recompute == True and os.path.exists(out_filename):
        os.remove(out_filename)
    if not os.path.exists(out_filename):
        align           = AlignIO.read(fasta, "fasta")
        align_columns   = []
        rec_names = []
        seq_index = 0
        for record in align:
            seq = str(record.seq)
            rec_names.append(record.id)
            index = 0
            for i in seq:
                if seq_index == 0:
                    align_columns.append(i)
                else:
                    align_columns[index] += i
                index += 1
                
            seq_index = 1

        align_clean = []
        seq_index = 0
        for col in align_columns:
            remove_col = False
            for x in lower_nucs:
                if x in col:
                    remove_col = True
                else:
                    continue
            if remove_col == False:
                index = 0
                for nuc in col:
                    if seq_index == 0:
                        align_clean.append(nuc)
                    else:
                        align_clean[index] += nuc
                    index += 1
            
            else:
                continue
            seq_index = 1
                    
        align1 = AlignIO.MultipleSeqAlignment([], Gapped(IUPAC.unambiguous_dna, "-"))
        for i, x in zip(align_clean, rec_names):
            temp_seq = SeqRecord(Seq(i, generic_dna), id = x)
            align1.append(temp_seq)
        AlignIO.write(align1, out_filename, "fasta")
        #record.seq = new_seq
        #AlignIO.write(align, out_filename, "fasta")
    return(out_filename)



def clean_to_phylip(align_outfile_name, phylip_filename, recompute):
        ## convert to phylip and remove redundant  sequences
   
    redundant_seqs_file = 'redundant_sequences.txt'
    if recompute == True and os.path.exists(phylip_filename):
        os.remove(phylip_filename)
        os.remove(redundant_seqs_file)
    if not os.path.exists(phylip_filename):
        try:
            red_seq = open(redundant_seqs_file, 'w')
            align           = AlignIO.read(align_outfile_name, "fasta")
            uniq_seqs = {} ## sequence is the key, values are the names
            ## make the alignment only have unique sequences
            for record in align:
                seq = str(record.seq)
                sid = str(record.id)
                ## remove empty records that didn't align
                if set(seq) == set("-"):
                    continue
                ## concatenate id's of records that are all matches
                if not seq in uniq_seqs:
                    uniq_seqs[seq] = sid
                else:
                    line_out = "{} {}".format(uniq_seqs[seq], sid)
                    red_seq.write(line_out)
                    continue
            red_seq.close()
            align1 = AlignIO.MultipleSeqAlignment([], Gapped(IUPAC.unambiguous_dna, "-"))
            for x in uniq_seqs:
                align1.add_sequence(uniq_seqs[x], x)
            if align1.get_alignment_length() >= 30:
                AlignIO.write(align1, phylip_filename, "phylip-sequential")
        except:
            print(align_outfile_name, 'too short to continue after trimming')
    return(phylip_filename)


def main_cds_tree(fasta, threads, recompute, id_map_filename = False, mafft = False):
    if threads == 'AUTO':
        tr = '1'
    else:
        tr = threads
    if mafft == True:
        align_outfile_name = run_mafft(fasta = fasta, threads = tr, recompute = recompute)
    else:
        fasta                = check_nucleotides_dialign(fasta)
        align_outfile_name   = run_dialign(fasta = fasta, recompute = recompute)
        align_outfile_name   = clean_dialign_fa(fasta = align_outfile_name, recompute = recompute)
    else:
            fasta              = clean_fasta_codons(fasta = fasta, recompute = recompute)
            align_outfile_name = run_prank(fasta = fasta, recompute = recompute)
        
    align_outfile_name   = run_trimal(afa = align_outfile_name, recompute = recompute)
    phylip_filename      = re.sub("afa.clean.trimmed", "phylip", align_outfile_name)
    phylip_filename      = clean_to_phylip(align_outfile_name = align_outfile_name, phylip_filename = phylip_filename, recompute = recompute)
            
    if '-yn' in sys.argv:
        phylip_yn00_filename = clean_phylip_yn00(align_outfile_name = align_outfile_name, recompute = recompute)
        run_yn00(phylip = phylip_yn00_filename, recompute = recompute, id_map = id_map_filename)
    if '-a' in sys.argv:
        sys.exit()
        
    ## align the phylip file
    if '-fasttree' in sys.argv:
        tree_filename = run_fasttree(fasta = align_outfile_name, threads = threads, recompute = recompute, cds = True)
    else:
        ## run modeltest to get the model
        output_file = 'modeltest/' + os.path.basename(align_outfile_name) + ".COMPLETE"
        if not recompute and os.path.exists(output_file):
            with open(output_file, 'r') as f:
                model = parse_model_test(in_file = output_file)
        else:
            tree_filename = run_iqtree(phylip = phylip_filename, threads = threads, recompute = recompute, cds = True)
    return(tree_filename)

    
def main_protein_tree(fasta, threads, recompute):
    mafft_outfile_name = run_mafft(fasta = fasta, threads = threads, recompute = recompute)
    align_outfile_name = run_trimal(afa = mafft_outfile_name, recompute = recompute)
    phylip_filename    = clean_to_phylip(align_outfile_name = align_outfile_name, recompute = recompute)
    #print(phylip_filename)
    if '-a' in sys.argv:
        sys.exit()
        
    ## align the phylip file
    if '-fasttree' in sys.argv:
        tree_filename = run_fasttree(phylip = phylip_filename, threads = threads, recompute = recompute)
    else:
        tree_filename = run_iqtree(phylip = phylip_filename, threads = threads, recompute = recompute, cds = False)
    return(tree_filename)
        


def main(argv):
    if len(argv) < 2:
        sys.stderr.write('''Usage: %s <fastafile> <threads> <-r>
        <fastafile>     the fasta file of orthologs to make a tree of
        <threads>       the number of threads to use
        <-r>            recompute everything
        <-a>            stop after aligning
        <-yn>           run yn00 to get dn/ds for cds sequences\n''' % os.path.basename(sys.argv[0]))
        sys.exit(1)
    fasta   = sys.argv[1]
    threads = sys.argv[2]
    if '-r' in sys.argv:
        recompute = True
    else:
        recompute = False
    if '-mafft' in sys.argv:
        mafft = True
    else:
        mafft = False
    if '-dialign' in sys.argv:
        dialign = True
    else:
        dialign = False
    print('Processing:', fasta)
        
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
    elif 'id_map.txt' in os.listdir('.'):
            id_map_filename = 'id_map.txt'
    else:
        id_map_filename = False

    tree_filename = main_protein_tree(fasta = fasta, threads = threads, recompute = recompute)
        
    ## pull the tree in and fix the newick back to names from the original fasta file
    if not os.path.exists('final_trees'):
        os.mkdir('final_trees')
        
    if id_map_filename:
        nwk_out = rename_newick(id_map = id_map_filename, nwk_in = tree_filename)
    else:
        nwk_out = os.path.join('final_trees', os.path.basename(tree_filename))
        shutil.copyfile(tree_filename, nwk_out)
    #print('Final tree:', nwk_out)

if __name__ == '__main__':
    main(sys.argv[1:])

