#!/usr/bin/env python

from BCBio import GFF
from Bio import SeqIO
import sys
import os, re

def main(argv):
    if len(argv) != 2 or '-h' in sys.argv:
        sys.stderr.write("Usage: %s <in.gbk> <genome_id>\n" % os.path.basename(sys.argv[0]))
        sys.exit(1)

    in_file = sys.argv[1]
    genome_id = sys.argv[2]

    in_handle = open(in_file)

    parent = ''
    gene   = ''
    for record in SeqIO.parse(in_handle, "genbank"):
        scaffold = record.id
        for feature in record.features:
            end = feature.location.end
            start = feature.location.start + 1
            strand = feature.location.strand
            if strand == -1:
                strand = '-'
            elif strand == 1:
                strand = '+'
            else:
                strand = '*'
            f_type = feature.type
            quals = feature.qualifiers
            #print(quals.keys())
            if 'score' in quals.keys():
                score = quals['score'][0]
            else:
                score = "-"
            if 'phase' in quals.keys():
                phase = quals['phase'][0]
            else:
                phase = '-'
            if f_type == 'cluster':
                parent = quals['note'][0]
                parent = 'cluster_' + parent.split(" ")[-1]
                attributes = "ID={0};Name={0};Note=contig_edge:{1},product:{2},GenomeID:{3}".format(parent, quals['contig_edge'][0], quals['product'][0], genome_id)
                source = 'antismash4'
            else:
                if f_type == "CDS":
                    if 'note' in quals.keys():
                        if 'gene' in quals.keys():
                            attributes = 'ID={0};Name={0};Parent={1};Note='.format(quals['gene'][0], parent)
                        elif 'locus_tag' in quals.keys():
                            attributes = 'ID={0};Name={0};Parent={1};Note='.format(quals['locus_tag'][0], parent)
                        else:
                            print('Unknown key, exiting', quals.keys())
                            sys.exit()
                    else:
                        if 'gene' in quals.keys():
                            attributes = 'ID={0};Name={0};Parent={1};Note='.format(quals['gene'][0], parent)
                        elif 'locus_tag' in quals.keys():
                            attributes = 'ID={0};Name={0};Parent={1};Note='.format(quals['locus_tag'][0], parent)
                        else:
                            print('Unknown key, exiting', quals.keys())
                            sys.exit()
                    if 'source' in quals.keys():
                        source = quals['source'][0]
                    else:
                        source = 'maker'
                    if 'gene' in quals.keys():
                        gene = quals['gene'][0]
                    else:
                        gene = quals['locus_tag'][0]
                elif f_type in ["CDS_motif", "PFAM_domain", "aSDomain"]:
                    source = quals['detection'][0]
                    gene = quals["locus_tag"][0]
                    if f_type == 'PFAM_domain':
                        attributes = "ID={0};Name={0};Parent={1};Dbxref={2};Note=description:{3}".format(quals['domain'][0], quals['locus_tag'][0], quals['db_xref'][0], quals['description'][0])
                    elif f_type == 'aSDomain':
                        attributes = "ID={0};Name={0};Parent={1};Note=".format(quals['domain'][0], quals['locus_tag'][0])
                    elif f_type == 'CDS_motif':
                        attributes = "ID={0};Name={0};Parent={1};Note=description:{2}".format(quals['motif'][0], quals['locus_tag'][0], quals['note'][0])
                else:
                    print("Missing feature, exiting program")
                    sys.exit()
                attributes += ",GenomeID:" + genome_id
                scaffold = re.sub(" .*$", "", scaffold)
                        
            items_all = [scaffold, source, f_type, str(start), str(end), str(score), str(strand), phase, attributes]
            print("\t".join(items_all))
    in_handle.close()

if __name__ == '__main__':
    main(sys.argv[1:])
