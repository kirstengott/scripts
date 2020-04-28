#!/usr/bin/env python3



import sys, os, re, pysam


## only keep orthogroups with dn/ds calculations
ortho_keep = [re.sub("\..*$", "", (x)) for x in os.listdir('yn_out') if 'full' not in x]


keep_cds = {}
out_dir = 'ani_cds_db'
if not os.path.exists(out_dir):
    os.mkdir(out_dir)


for x in os.listdir('ortho_fastas'):
    name = re.sub('\..*$', '', x)
    path = os.path.join('ortho_fastas', x)
    if name in ortho_keep:
        with pysam.FastxFile(path) as fh:
            for entry in fh:
                if entry.name in keep_cds:
                    keep_cds[entry.name].append(entry.sequence)
                else:
                    keep_cds[entry.name] = [entry.sequence]
    else:
        continue

outstring = ">{}\n{}\n"
for x in keep_cds:
    file_out = os.path.join(out_dir, x)
    file_out = file_out + '.fa'
    seq_count = 1
    with open(file_out, 'w') as fh:
        for seq in keep_cds[x]:
            seq_name = 'seq' + str(seq_count)
            outline = outstring.format(seq_name, seq)
            fh.write(outline)
            seq_count += 1
    seq_count = 0
            
        
        
    
        
                
