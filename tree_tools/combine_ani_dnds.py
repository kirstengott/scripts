#!/usr/bin/env python3
import sys, os, re, pysam




dn_ds = open('dn_ds_summary.txt', 'r')
ani = open('fastani/ANIvis.tsv', 'r')

vals = {}

na_count = 0
total_count = 0
for line in dn_ds:
    line = line.strip().split(',')
    g1 = re.sub("\..*$", "", line[0])
    g2 = re.sub("\..*$", "", line[1])
    if g1 == g2:
        continue
    c1 = g1 + "," + g2
    val = line[2]
    vals[c1] = {}
    vals[c1]['dn_ds'] = val

header = 0
for line in ani:
    if header == 0:
        header = 1
        continue
    line = line.strip().split()
    g1 = re.sub("\..*$", "", os.path.basename(line[0]))
    g2 = re.sub("\..*$", "", os.path.basename(line[1]))
    if g1 == g2:
        continue
    c1 = g1 + "," + g2
    c2 = g2 + "," + g1
    an = float(line[2])
    in_dict = 0
    if c1 in vals:
        in_dict = 1
    elif c2 in vals:
        in_dict = 2
    if in_dict == 1:
        if 'ani' in vals[c1]:
            vals[c1]['ani'] =  (vals[c1]['ani'] + an)/2
        else:
            vals[c1]['ani'] = an
    elif in_dict == 2:
        if 'ani' in vals[c2]:
            vals[c2]['ani'] =  (vals[c2]['ani'] + an)/2
        else:
            vals[c2]['ani'] = an

dn_ds.close()
ani.close()

with open("ani_dn_ds_summary.txt", 'w') as fh:
    for i in vals:
        string = "{},{},{}\n".format(i, vals[i]['dn_ds'], str(vals[i]['ani']))
        fh.write(string)
