#!/usr/bin/env python3



import sys, os, re, pysam

def main(argv):
    if len(argv) <2:
        sys.stderr.write('''Usage: %s <yn_out> <fasta_dir>
        <yn_out>       path to yn_out directory
        <fasta_dir>    path to input fasta file directory \n''' % os.path.basename(sys.argv[0]))
        sys.exit()

    yn_out = sys.argv[1]
    fa     = sys.argv[2]

    ## only keep orthogroups with dn/ds calculations
    ortho_keep = [re.sub("\..*$", "", (x)) for x in os.listdir(yn_out) if 'full' not in x]


    keep_cds = {}
    out_dir = 'ani_cds_db'
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)


    for x in os.listdir(fa):
        name = re.sub('\..*$', '', x)
        path = os.path.join(fa, x)
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
            
        
        
    
        
                
if __name__ == '__main__':
    main(sys.argv[1:])
