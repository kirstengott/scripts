import sys, os, re

def rename_newick(id_map, nwk_in):
    print(nwk_in)
    nwk_f = open(nwk_in, 'r')
    id_map_f = open(id_map, 'r')
    nwk = nwk_f.readlines()[0]
    nwk_f.close()
    for line in id_map_f:
        line = line.rstrip().split()
        temp_id = line[0] + ':'
        new_id = re.sub(",", ".", line[1]) + ":"
        nwk = re.sub(temp_id, new_id, nwk)
    id_map_f.close()
    if not os.path.exists('final_trees'):
        os.mkdir('final_trees')
    nwk_out = os.path.join('final_trees', os.path.basename(nwk_in))
    nwk_out_h = open(nwk_out, 'w')
    nwk_out_h.write(nwk)
    nwk_out_h.close()
    return(nwk_out)


rename_newick(nwk_in = sys.argv[1], id_map = sys.argv[2])
