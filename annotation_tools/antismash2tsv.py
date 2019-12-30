#!/usr/bin/env python

from BCBio import GFF
from Bio import SeqIO
import argparse


parser = argparse.ArgumentParser(description="")
parser.add_argument("-i", help = "path to input genbank")
parser.add_argument('-t', help = 'path to geneclusters/txt')
args = parser.parse_args()

in_file = args.i


in_handle = open(in_file)


txt_file = args.t

gene_names_dict = {}
chr_names_dict = {}


for record in SeqIO.parse(in_handle, "genbank"):
    chr_names_dict[record.id] = record.description
    for feature in record.features:
        if "translation" in feature.qualifiers and feature.type == "CDS":
            locus_tag = feature.qualifiers["locus_tag"][0]
            mrna_id = feature.qualifiers["Parent"][0]
            gene_names_dict[locus_tag] = [mrna_id]
               
lines = [x.strip().split("\t") for x in open(txt_file, 'r')]


sub_lines = [[x[0], x[2], x[3]] for x in lines]


#print("\t".join(['chromosome', 'metabolite', 'mRNA']))

for l in sub_lines:
    chr = chr_names_dict[l[0]]
    metabolite = l[1]
    gene_tags = l[2].split(";")
    #cds = []
    genes = []
    for gene in gene_tags:
        locus_info = gene_names_dict[gene]
        #cds.append(locus_info[0])
        genes.append(locus_info[0].strip())
        #cds_all = ",".join(set(cds))
    genes_all = ",".join(genes)
    print("\t".join([chr, metabolite, genes_all]))
    
    
in_handle.close()

