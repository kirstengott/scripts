#!/usr/bin/env python

from BCBio import GFF
from Bio import SeqIO
import sys
import os

def main(argv):
    if len(argv) != 3:
        sys.stderr.write("Usage: %s <in.gbk> <in.geneclusters.txt> <genome_id>\n" % os.path.basename(sys.argv[0]))
        sys.exit(1)

    in_file = sys.argv[1]
    txt_file = sys.argv[2]
    genome_id = sys.argv[3]

    in_handle = open(in_file)

    scaffold_names_dict = {}
    
    with open(txt_file, 'r') as t_file:
        for line in t_file:
            l = line.strip().split("\t")
            all_antismash_genes = l[3].split(";")
            for gene in all_antismash_genes:
                scaffold_names_dict[gene] = [l[1], l[2]]
    gene = ''
    gff3_lines = []
    for record in SeqIO.parse(in_handle, "genbank"):
        for feature in record.features:
            end = feature.location.end
            start = feature.location.start
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
                cluster_items = [source, f_type, str(start), str(end), str(score), str(strand), phase, attributes]
            else:
                if f_type == "CDS":
                    if 'note' in quals.keys():
                        attributes = 'ID={0};Name={0};Parent={1};Note='.format(quals['gene'][0], parent)
                    else:
                        attributes = 'ID={0};Name={0};Parent={1};Note='.format(quals['gene'][0], parent)
                    source = quals['source'][0]
                    gene = quals['gene'][0]
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
                if gene:
                    if gene in scaffold_names_dict.keys():
                        scaffold = scaffold_names_dict[gene][0]
                    else:
                        scaffold = 'unk'
                attributes += ",GenomeID:" + genome_id
                items_all = [scaffold, source, f_type, str(start), str(end), str(score), str(strand), phase, attributes]
                gff3_lines.append(items_all)
                gene = ''

    cluster_items.insert(0,scaffold)
    gff3_lines.insert(0, cluster_items)

    all_lines = "\n".join(["\t".join(x) for x in gff3_lines])

    print(all_lines)

    in_handle.close()

if __name__ == '__main__':
    main(sys.argv[1:])
