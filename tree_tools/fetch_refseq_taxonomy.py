#!/usr/bin/env python3

from subprocess import check_output
import sys
import os
import shutil
from gbk2faa import gbk2faa
import re
#from pathlib import Path


## make this stat the file to see when it was pulled down
# if os.path.exists('~/refseq'):
#     os.remove('~/refseq')
#     os.remove(ncbi)

    
#os.system('ncbi-genome-download -F protein-fasta -l all -p 90 -m ncbi_metatdata.txt -v fungi')
#os.system('ncbi-genome-download -s genbank -F genbank,protein-fasta -l all -p 90 -m genbank_metatdata.txt -v fungi')


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
    except:
        print('no proteins in file')
    os.system('gzip ' + gbk)
    if os.stat(gbk_faa).st_size == 0:
        return(False)
        os.remove(gbk_faa)
    else:
        return(gbk_faa)

def get_taxonomy(ncbi, taxa_pull, output_directory):
    taxa_pull = [x.lower() for x in taxa_pull]
    ncbi_f = open(ncbi, 'r')
    keep = ['assembly_accession', 'taxid', 'organism_name', 'local_filename', 'assembly_level', 'gbrs_paired_asm']
    if not os.path.exists(output_directory):
        os.mkdir(output_directory)
    outfile = open(os.path.join(output_directory, 'taxonomy.txt'), 'a')
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
            print(genome_info['organism_name'])
            command = 'efetch -db taxonomy -id '+ taxid +' -format xml | xtract -pattern Taxon -block "*/Taxon" -tab "\t" -element ScientificName'
            out = check_output(command, shell = True)
            taxonomy = out.decode("utf-8").strip().split()
            taxa_lower = [x.lower() for x in taxonomy]
            check_taxa = list(set([True for taxa in taxa_pull if taxa in taxa_lower]))
            ## if there are no desired taxa in this entry, continue to the next iteration
            if len(check_taxa) > 0:
                check_taxa = check_taxa[0]
                #print(genome_info['organism_name'])
            else:
                continue
            if check_taxa:
                if genome_info['local_filename'].endswith('faa.gz'):
                    ## only keeping genomes that are annotated
                    outfile_name = genome_info['assembly_accession'] + ".faa.gz" ## rename as the genbank name
                    shutil.copyfile(genome_info['local_filename'], os.path.join(output_directory, outfile_name))
                    genome_info['local_filename'] = outfile_name
                # elif genome_info['local_filename'].endswith('gbff.gz'):
                #     ## what to do if the file is from genbank
                #     ## create a consistent naming schema, if there is a refseq name, use it, otherwise use the genbank name.
                #     if genome_info['gbrs_paired_asm'] == 'na':
                #         outfile_name = genome_info['assembly_accession']
                #     else:
                #         outfile_name = genome_info['gbrs_paired_asm']
                        
                #     direct = os.path.dirname(genome_info['local_filename'])
                #     proteins_annotated = [x for x in os.listdir(direct) if "_protein.faa.gz" in x]
                #     if len(proteins_annotated) == 0:
                #         local_proteins = check_gbk_faa(genome_info['local_filename'])
                #         if local_proteins:
                #             outfile_name = outfile_name + ".faa"
                #             print(outfile_name)
                #             shutil.copyfile(local_proteins, os.path.join(output_directory, outfile_name))
                #     else:
                #         proteins_annotated = proteins_annotated[0]
                #         outfile_name = outfile_name + ".faa.gz"
                #         print(outfile_name)
                #         shutil.copyfile(os.path.join(direct, proteins_annotated), os.path.join(output_directory, outfile_name))

                genome_info['local_filename'] = outfile_name
                out = "{},{}\n".format(",".join([genome_info[x] for x in keep]), ";".join(taxonomy))
                outfile.write(out)


    outfile.close()
    ncbi_f.close()


def usage():
    sys.stderr.write("Usage: %s <output_dir> <metadata_file> <taxa_grab>\n" % os.path.basename(sys.argv[0]))
    sys.exit(1)
    
def main(argv):
    if len(argv) < 2:
        usage()
    elif "-h" in sys.argv[1]:
        usage()
    output_directory = sys.argv[1]
    metadata = sys.argv[2]
    if len(argv) == 3:
        taxa_pull = ['hypocreales', 'Aspergillus', 'Saccharomyces'] ## only pull representatives for outgroups to hypocreales
    else:
        taxa_pull = ['hypocreales']
    get_taxonomy(ncbi = metadata, taxa_pull = taxa_pull, output_directory = output_directory)
    
if __name__ == "__main__":
    main(sys.argv[1:])
