#!/usr/bin/env python3
import sys


def main(argv):
    if len(argv) != 2:
        sys.stderr.write("Usage: %s <targetp_outfile> <id_map>\n" % os.path.basename(sys.argv[0]))
        sys.exit(1)
    f        = open(sys.argv[1], 'r')
    ids_file = open(sys.argv[2], 'r')

    ids_diction = {}
    for line in ids_file:
        l = line.strip().split("\t")
        ids_diction[l[0]] = l[1]

    #print(ids_diction)
    count = 0 
    for line in f:
        liner = line.strip().split()
        if 'cutoff' in line:
            pass
        elif '--' in line:
            pass
        elif count > 7:
            new_id = ids_diction[liner[0]]
            new_line = new_id + "\t" + "\t".join(liner[1:])
            print(new_line)
        count += 1
    ids_file.close()
    f.close()

if __name__ == '__main__':
    main(sys.argv[1:])
