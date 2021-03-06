from Bio import AlignIO

import sys

align = AlignIO.read(sys.argv[1], "phylip")

out = sys.argv[2]

with open(out, 'w') as fh:
    fh.write("seq gaps matches proportion")


    for i in align:
        seq = i.seq
        gaps = seq.count("-")
        matches = seq.count('a') + seq.count('t') + seq.count('g') + seq.count('c')
        prop = (gaps/(matches+gaps))*100
        fh.write("{} {} {} {}".format(i.id, gaps, matches, prop))
        if prop > 30:
            print(i.id, prop)
        
