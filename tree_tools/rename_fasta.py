import sys, os, re

def rename_fasta(id_map, fa):
    id_map_f = open(id_map, 'r')
    id_map = {}
    for line in id_map_f:
        line = line.rstrip().split()
        id_map[line[0]] = line[1]

    id_map_f.close()
    
    nwk_out = os.path.basename(fa) + ".named"
    nwk_out_h = open(nwk_out, 'w')

    f = open(fa, 'r')
    for line in f:
        if line.startswith(">"):
            line = line.strip()
            t_id = re.sub(">", "", line)
            new_id = ">" + id_map[t_id] + "\n"
            nwk_out_h.write(new_id)
        else:
            nwk_out_h.write(line)

    nwk_out_h.close()
    return(nwk_out)


rename_fasta(fa = sys.argv[1], id_map = sys.argv[2])
