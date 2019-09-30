#!/usr/bin/python

import re
import os.path
import argparse

parser = argparse.ArgumentParser(description='Fix the maker GFF to make it easier for antismash.')
parser.add_argument('gff',
                    type=str,
                    nargs='+',
                    help='the gff to fix')

args = parser.parse_args()

file = args.gff[0]
file_base = os.path.splitext(file)[0]


new_gff  = file_base + "_modified.gff"
map_file = file_base + "_id_map.txt"



n_gff  = open(new_gff, 'w')
m_file = open(map_file, 'w')
f      = open(file, 'r')




id_map = []

desc_ini = ''
id_counter = 0

for line in f:
    liner = str(line)
    items = line.strip().split("\t")
    ## define logic to remove comment lines
    if len(items) > 1:
        ## pull out the id section of the description field
        desc = items[8].split(";")[0].strip()
        desc_sub = re.sub('ID=', '', desc)
        ## define logic to group gene id records and initiate ID renaming for the first gene
        if desc_ini == '':
            desc_ini = desc_sub
        if desc_ini in desc_sub:
            new_id = 'id' + str(id_counter)
            new_line = re.sub(pattern=desc_ini, repl=new_id, string=liner)
            n_gff.write(new_line)
            id_mapper = desc_ini + "\t" + new_id + "\n"
            if id_mapper not in id_map:
                id_map.append(id_mapper)
        else:
            desc_ini = desc_sub
            id_counter += 1
    else:
        pass

## write out the key to the gff ids
[m_file.write(x) for x in id_map]    


## close the files
f.close()
m_file.close()
n_gff.close()
