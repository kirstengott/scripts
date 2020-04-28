#!/usr/bin/env python3

import sys, os, re, pysam


## only keep orthogroups with dn/ds calculations
ortho_keep = [os.path.join('yn_out', x) for x in os.listdir('yn_out') if 'full' not in x]
outfile = 'dn_ds_summary.txt'
ds_vals = {}

na_count = 0
total_count = 0
for x in ortho_keep:
    #name = os.path.splitext(re.sub('\..*$', '', os.path.basename(x)))
    header = 0
    with open(x, 'r') as fh:
        for line in fh:
            if header == 0:
                header = 1
                continue
            else:
                line = line.strip().split(',')
                g1 = os.path.splitext(line[0])[0]
                g2 = os.path.splitext(line[1])[0]
                c1 = g1 + "," + g2
                c2 = g2 + "," + g1
                try:
                    dn_ds = float(line[2])/float(line[4])
                    total_count += 1
                except:
                    dn_ds = 'NA'
                    na_count +=1
                in_dict = 0
                if c1 in ds_vals:
                    in_dict = 1
                elif c2 in ds_vals:
                    in_dict = 2
                else:
                    ## put c1 in first
                    ds_vals[c1] = [dn_ds]
                    
                if in_dict == 1:
                    ds_vals[c1].append(dn_ds)
                elif in_dict == 2:
                    ds_vals[c2].append(dn_ds)

of = open(outfile, 'w')
for i in ds_vals:
    na_remove = [x for x in ds_vals[i] if x != 'NA']
    if len(na_remove) > 0:
        dn_ds = sum(na_remove)/len(na_remove)
        outline = i + "," + str(dn_ds) + "\n"
        of.write(outline)
    else:
        num_na += 1
        
of.close()
