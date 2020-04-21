#!/usr/bin/env python3
from subprocess import check_output
import sys
import os
import shutil
from gbk2faa import gbk2faa
import re




def run_ncbi_genome_download(threads, ncbi_type):
    if ncbi_type == 'refseq':
        refseq_command = 'ncbi-genome-download -F protein-fasta -l all -p {} -m refseq_metadata.txt -v fungi'.format(threads)
        os.system(refseq_command)
    elif ncbi_type == 'genbank':
        genbank_command = 'ncbi-genome-download -s genbank -F protein-fasta -l all -p {} -m genbank_metadata.txt -v fungi'.format(threads)
        os.system(genbank_command)
    else:
        print('unknown recompute type, exiting')
        sys.exit(1)

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

def parse_metadata(ncbi):
    ''' Input: the metadata table from ncbi-genome-download
        Output: a list of dictionaries of all of the lines
    '''
    ncbi_f = open(ncbi, 'r')
    keep = ['assembly_accession', 'taxid', 'organism_name', 'local_filename', 'assembly_level', 'gbrs_paired_asm']
    all_lines = []
    header = 0 
    for line in ncbi_f:
        l = line.rstrip().split("\t")
        if header == 0:
            keep_cols = [l.index(x) for x in keep]
            header = 1
        else:
            genome_items = [l[x] for x in keep_cols]
            genome_info = dict(zip(keep, genome_items))
            all_lines.append(genome_info)
    ncbi_f.close()
    return(all_lines)

def efetch_taxa(taxid):
    ''' Input: comma seperated list of ncbi taxids
        Output: a dictionary of lists of taxonomies for each taxid
        Method: Runs the ncbi 'efetch' command on the list of taxids.
    '''
    command = 'efetch -db taxonomy -id '+ taxid +' -format xml | xtract -pattern Taxon -block "*/Taxon" -tab "\t" -element ScientificName'
    out = check_output(command, shell = True)
    taxonomy = out.decode("utf-8").split("\n")
    taxonomy = [x.lower().split("\t") for x in taxonomy]
    taxa_list = taxid.split(",")
    taxonomy_dict = dict(zip(taxa_list, taxonomy))
    return(taxonomy_dict)


def get_taxonomy(ncbi, taxa_pull, output_directory):
    taxa_pull = [x.lower() for x in taxa_pull]
    if not os.path.exists(output_directory):
            os.mkdir(output_directory)
    outtable_name = os.path.join(output_directory, 'taxonomy.txt')
    outfile = open(outtable_name, 'a')                        
    all_lines = parse_metadata(ncbi)
    all_taxa_ids = ",".join([x['taxid'] for x in all_lines])
    all_taxa = efetch_taxa(all_taxa_ids)
    #keep = all_lines[0].keys()
    for i in all_lines:
        genome_info = i
        taxid = genome_info['taxid']
        taxonomy = all_taxa[taxid]
        check_taxa = list(set([True for taxa in taxa_pull if taxa in taxonomy]))
        ## if there are no desired taxa in this entry, continue to the next iteration
        if len(check_taxa) > 0:
            check_taxa = check_taxa[0] ## Flag for whether or not we want this taxon
        else:
            continue
        if check_taxa:
            if 'faa' in genome_info['local_filename']:
                ext = '.faa.gz'
            elif 'cds_from_genomic' in genome_info['local_filename']:
                ext = '.cds_from_genomic.fna.gz'
            else:
                print(genome_info['local_filename'], "not included, fix script")
                sys.exit()
            ## only keeping genomes that are annotated
            if 'genbank' in ncbi:
                outfile_name = genome_info['assembly_accession'] + ext ## rename as the genbank name
                keep = ['assembly_accession', 'taxid', 'organism_name', 'local_filename', 'assembly_level', 'gbrs_paired_asm']

            elif 'refseq' in ncbi:
                outfile_name = genome_info['gbrs_paired_asm'] + ext ## rename as the genbank name
                keep = ['gbrs_paired_asm', 'taxid', 'organism_name', 'local_filename', 'assembly_level', 'assembly_accession']
            else:
                print("Couldn't detect database type from metadata filename, using default outputfile behavior")
                outfile_name = genome_info['assembly_accession'] + ext ## rename as the genbank name
                keep = ['assembly_accession', 'taxid', 'organism_name', 'local_filename', 'assembly_level', 'gbrs_paired_asm']
                

            shutil.copyfile(genome_info['local_filename'], os.path.join(output_directory, outfile_name))
            genome_info['local_filename'] = outfile_name
            out = "{},{}\n".format(",".join([genome_info[x] for x in keep]), ";".join(taxonomy))
            outfile.write(out)
    outfile.close()

    ## cleanup the output file
    infile = open(outtable_name, 'r')
    outtable2_name = os.path.join(output_directory, 'taxonomy_temp.txt')
    outfile = open(outtable2_name, 'w')
    unique_lines = []
    for line in infile:
        if line not in unique_lines:
            unique_lines.append(line)
            outfile.write(line)
        else:
            continue
    outfile.close()
    infile.close()
    os.remove(outtable_name)
    shutil.move(outtable2_name, outtable_name)




def usage():
    usage = '''Usage: %s <output_dir> <metadata_file> <recompute> <taxa_grab>
               <output_dir> where to house the fasta files and taxonomy table
               <metadata_file> which metadata to use
               <recompute> either False or the number of threads to redownload with
               <taxa_grab> if supplied will pull down supplied taxa\n''' % os.path.basename(sys.argv[0])
    sys.stderr.write(usage)
    sys.exit(1)
    
def main(argv):
    if len(argv) < 3:
        usage()
    elif "-h" in sys.argv[1]:
        usage()
    output_directory = sys.argv[1]
    metadata = sys.argv[2]
    recompute = sys.argv[3]
    ## fix this so it only removes the match to the metadata
    if recompute != 'False':
        which_rec = metadata.split("_")[0]
        outpath = os.path.join("~" + which_rec)
        if os.path.exists(outpath):
            os.remove(outpath)
            os.remove(metadata)
        run_ncbi_genome_download(threads = recompute, ncbi_type = which_rec)
    if len(argv) == 4:
        taxa_pull = [sys.argv[4]] ## only pull representatives for outgroups to hypocreales
    else:
        taxa_pull = ['hypocreales']
    get_taxonomy(ncbi = metadata, taxa_pull = taxa_pull, output_directory = output_directory)
    
if __name__ == "__main__":
    main(sys.argv[1:])
