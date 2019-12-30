#!/usr/bin/env python3

import sys
import os

def main(argv):
    if len(argv) != 2:
        sys.stderr.write("Usage: %s <interpro_all> <interpro_reduced>\n" % os.path.basename(sys.argv[0]))
        sys.exit(1)
    
    f = sys.argv[1]
    o = sys.argv[2]

    f_h = open(f, 'r')

    gene_info = {}

    all_tools = []

    for line in f_h:
        l = line.strip().split("\t")
        len_out = len(l) + 1
        ## grab interpro, go, and pathway terms if they are there


        ## add in logic for annotations that are only sometimes in the line
        if len(l) > 10:
            extra_annots = {}
            try:
                extra_annots['interpro'] = l[11]
                if 'interpro' not in all_tools:
                    all_tools.append('interpro')
            except:
                pass
            try:
                extra_annots['GO'] = l[13].split("|")
                #print(l[0])
                #print(extra_annots['GO'])
                if 'GO' not in all_tools:
                    all_tools.append('GO')
            except:
                pass
            try:
                extra_annots['PATHWAY'] = l[14]
                if 'PATHWAY' not in all_tools:
                    all_tools.append('PATHWAY')
            except:
                pass

            
        ## add unique tools to the list of tools
        if l[3] not in all_tools:
            all_tools.append(l[3])            
        if l[0] not in gene_info: ## add and initialize gene keys
            gene_info[l[0]] = {l[3] : [l[4]]}
        else:## append
            if l[3] not in gene_info[l[0]]:
                gene_info[l[0]][l[3]] = [l[4]]
            else:
                gene_info[l[0]][l[3]].append(l[4])

        ## add in extra annotations
        if extra_annots:
            for key in extra_annots:
                if key not in all_tools:
                    all_tools.append(key)

                ## all genes should be initialized at this point
                if key not in gene_info[l[0]]:
                    if key == 'GO':
                        gene_info[l[0]][key] = extra_annots[key]
                    else:
                        gene_info[l[0]][key] = [extra_annots[key]]
                else:
                    if key == 'GO':
                        for x in extra_annots[key]:
                            if x not in gene_info[l[0]][key]:
                                #print(gene_info[l[0]][key])
                                gene_info[l[0]][key].append(x)
                            else:
                                pass
                    #gene_info[l[0]][key].append(extra_annots[key])
    #print(all_tools)
    #print(gene_info['maker-NODE_141_length_61473_cov_48.100349-snap-gene-0.3-mRNA-1'])
    #sys.exit()       

                

    simplified_f = open(o, 'w')
    simplified_f.write('Gene' + "\t" + "\t".join(all_tools) + "\n")
    for gene in gene_info:
        sub_dict = gene_info[gene]
        annot = []
        for tool in all_tools:
            if tool in sub_dict:
                all_domains = []
                for i in sub_dict[tool]:
                    if i:
                        if i not in all_domains:
                            all_domains.append(i)
                out = ";".join(all_domains)
            else:
                out = 'NA'
            annot.append(out)
        simplified_f.write(gene + "\t" + "\t".join(annot) + "\n")
        

if __name__ == '__main__':
    main(sys.argv[1:])
