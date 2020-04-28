import os, sys, re, pysam, shutil

of_results = 'Results_Apr20'
seq_input = 'cds'
outdir = 'ortho_fastas'
seq_type = 'cds'

if not os.path.exists(outdir):
    os.mkdir(outdir)
else:
    shutil.rmtree(outdir)
    os.mkdir(outdir)

seqs = os.path.join(of_results, 'Orthogroup_Sequences')


gene_dict = {}

for i in os.listdir(seq_input):
    infile = os.path.join(seq_input, i)
    genome = os.path.splitext(i)[0]
    with pysam.FastxFile(infile) as fh:
        for entry in fh:
            if seq_type == 'cds':
                if entry.comment:
                    if 'protein_id' in entry.comment:
                        name = [x for x in entry.comment.split() if 'protein_id' in x][0]
                        name = re.sub('\[protein_id=', '', name)
                        name = re.sub('\]', '', name)
                        gene_dict[name] = [genome, entry.sequence]
                    else:
                        ## these are likely pseudo genes from my tests, so we will omit them
                        continue
                else:
                    gene_dict[entry.name] = [genome, entry.sequence]
            else:
                gene_dict[entry.name] = genome

out_format = ">{}\n{}\n"
for i in os.listdir(seqs):
    infile = os.path.join(seqs, i)
    outfile = os.path.join(outdir, i)
    out_f = open(outfile, 'w')
    unique_genomes = set()
    with pysam.FastxFile(infile) as fh:
        for entry in fh:
            if seq_type == 'cds':
                name = gene_dict[entry.name][0]
                seq  = gene_dict[entry.name][1]
                outstring = out_format.format(name, seq)
            else:
                name = gene_dict[entry.name] 
                outstring = out_format.format(name, entry.sequence)
            unique_genomes.add(name)
            out_f.write(outstring)
    out_f.close()
    if len(unique_genomes) < 3:
        os.remove(outfile)
    












