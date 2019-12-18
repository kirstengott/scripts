#!/usr/bin/env python3

from subprocess import check_output
import sys
import os
import shutil
from gbk2faa import gbk2faa
import re

ncbi = 'ncbi_metadata.txt'
genbank = 'genbank_metadata.txt'

## make this stat the file to see when it was pulled down
# if os.path.exists('~/refseq'):
#     os.remove('~/refseq')
#     os.remove(ncbi)

    
#os.system('ncbi-genome-download -F protein-fasta -l all -p 90 -m ncbi_metatdata.txt -v fungi')
#os.system('ncbi-genome-download -s genbank -F genbank -l all -p 90 -m genbank_metatdata.txt -v fungi')



taxa_pull_refseq = ['hypocreales', 'Aspergillus', 'Saccharomyces'] ## only pull representatives for outgroups to hypocreales
taxa_pull_genbank = ['hypocreales']



output_directory = 'hypocreales_tree'

def check_gbk_faa(gbk):
    ## if the genbank file has CDS sequences, this function pulls
    ## out the proteins and returns the file name of the proteins file
    ## otherwise returns 'False'
    gbk_dir = os.path.dirname(gbk)
    os.system('gunzip ' + gbk)
    gbk_r = re.sub('.gz', '', gbk)
    gbk_faa = re.sub('gbff', 'faa', gbk_r)
    try:
        gbk2faa(gbk_filename = gbk_r, faa_filename = gbk_faa)
        return(gbk_faa)
    except:
        return(False)
    # if os.stat(gbk_faa).st_size == 0:
    #     return(False)
    # else:
    #     return(gbk_faa)

def get_taxonomy(ncbi, taxa_pull, output_directory):
    taxa_pull = [x.lower() for x in taxa_pull]
    ncbi_f = open(ncbi, 'r')
    keep = ['assembly_accession', 'taxid', 'organism_name', 'local_filename', 'assembly_level', 'gbrs_paired_asm']
    if not os.path.exists(output_directory):
        os.mkdir(output_directory)
    outfile = open(os.path.join(output_directory, 'ncbi_taxonomy.txt'), 'a')
    header = 0 
    for line in ncbi_f:
        l = line.rstrip().split("\t")
        if header == 0:
            keep_cols = [l.index(x) for x in keep]
            header = 1
        else:
            genome_items = [l[x] for x in keep_cols]
            genome_info = dict(zip(keep, genome_items))

            taxid = genome_info['taxid']
            command = 'efetch -db taxonomy -id '+ taxid +' -format xml | xtract -pattern Taxon -block "*/Taxon" -tab "\t" -element ScientificName'
            out = check_output(command, shell = True)
            taxonomy = out.decode("utf-8").strip().split()
            taxa_lower = [x.lower() for x in taxonomy]
            check_taxa = list(set([True for taxa in taxa_pull if taxa in taxa_lower]))
            ## if there are no desired taxa in this entry, continue to the next iteration
            if len(check_taxa) > 0:
                check_taxa = check_taxa[0]
                print(genome_info['organism_name'])
            else:
                continue
            if check_taxa:
                if genome_info['local_filename'].endswith('faa.gz'):
                    ## what to do if the file is from refseq
                    outfile_name = genome_info['gbrs_paired_asm'] + ".faa.gz"
                    shutil.copyfile(genome_info['local_filename'], os.path.join(output_directory, outfile_name))
                    genome_info['local_filename'] = outfile_name

                elif genome_info['local_filename'].endswith('gbff.gz'):
                    ## what to do if the file is from genbank
                    local_proteins = check_gbk_faa(genome_info['local_filename'])
                    if local_proteins:
                        outfile_name = genome_info['assembly_accession'] + ".faa"
                        print(outfile_name)
                        shutil.copyfile(local_proteins, os.path.join(output_directory, outfile_name))
                    
                genome_info['local_filename'] = outfile_name
                out = "{},{}\n".format(",".join([genome_info[x] for x in keep]), ";".join(taxonomy))
                outfile.write(out)


    outfile.close()
    ncbi_f.close()

    
#get_taxonomy(ncbi = ncbi, taxa_pull = taxa_pull_refseq, output_directory = output_directory)
get_taxonomy(ncbi = genbank, taxa_pull = taxa_pull_genbank, output_directory = output_directory)        
#efetch -db taxonomy -id 4754 -format xml | xtract -pattern Taxon -block "*/Taxon" -tab "\t" -element ScientificName
