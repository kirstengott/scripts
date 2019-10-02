#!/usr/bin/env python

import sys



def main(argv):
    if len(argv) != 2:
        sys.stderr.write("Usage: %s <wolfpsort.out> <full or partial> >output\n" % os.path.basename(sys.argv[0]))
        sys.exit(1)
    inp = sys.argv[1]
    for line in open(inp, 'r'):
        if '#' in line:
            pass
        else:
            ln = line.strip().split(",")
            gene = ln[0].split(" ")[0]
            first_item = ln[0].split(" ")[1]
            first_val = ln[0].split(" ")[2]
            if sys.argv[2] == 'full':
                print("{0}\t{1}\t{2}".format(gene, first_item, first_val))
                for i in ln[1:len(ln)]:
                    item = i.strip().split(" ")[0]
                    val = i.strip().split(" ")[1]
                    print("{0}\t{1}\t{2}".format(gene, item, val))
            else:
                print("{0}\t{1}".format(gene, first_item))


        
        
if __name__ == '__main__':
    main(sys.argv[1:])
