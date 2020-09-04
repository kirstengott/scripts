#!/usr/bin/env python3

import sys, os, re, pysam

def main(argv):
    if len(argv) <2:
        sys.stderr.write('''Usage: %s <yn_out> <out_prefix>
        <yn_out>       path to yn_out directory
        <out_prefix>   prefix to attach to output\n''' % os.path.basename(sys.argv[0]))
        sys.exit()

    yn_out = sys.argv[1]
    out_p  = sys.argv[2]
    
    ## only keep orthogroups with dn/ds calculations
    ortho_keep = [os.path.join(yn_out, x) for x in os.listdir(yn_out) if 'full' not in x]
    outfile = out_p + "-" + 'dn_ds_summary.txt'
    outfile_all = open(out_p + "-" + 'dn_ds_all.txt', 'w')
    ds_vals = {}

    na_count = 0
    total_count = 0
    for x in ortho_keep:
        name = os.path.splitext(re.sub('\..*$', '', os.path.basename(x)))[0]
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
                    line_all_out = "{},{},{},{},{},{}\n".format(name, g1, g2, dn_ds, line[2], line[4])
                    outfile_all.write(line_all_out)
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

    outfile_all.close()
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

if __name__ == '__main__':
    main(sys.argv[1:])
    
