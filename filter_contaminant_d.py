import pysam
import sys
import os
import re

file_in = sys.argv[1]
fasta = sys.argv[2]


t_in = 0
e_in = 0
d_in = 0

c_dict = {}
exclude = []
          
with open(file_in, 'r') as fh:
    for line in fh:
        line = line.strip()

        ### Starting logic to create contaminants to exclude
        if 'Exclude' in line:
                e_in = 1
                continue
        elif e_in == 1:
            e_in +=1 
            continue
        elif e_in == 2:
            lin = line.split("\t")
            if len(lin) < 3:
                e_in = 0
                continue
            exclude.append(lin[0])
            

        ### Starting logic to create dictionary of sequences to trim    
        elif 'Trim:' in line:
            t_in = 1
            continue
        elif t_in == 1:
            t_in += 1
        elif t_in == 2:
            ele = line.split("\t")
            if len(ele) < 4:
                t_in = 0
                continue
            if ele[0] not in c_dict:
                if "," in ele[2]:
                    item = ele[2].split(",")
                    items = [x.split("..") for x in item]
                else:
                    items = [ele[2].split("..")]
            else:
                continue
            c_dict[ele[0]] = items ## put the range in the dictionary

        ### Starting logic to create list of dupes to exclude
        
        elif 'Duplicated' in line:
            d_in = 1
            continue
        elif d_in == 1:
            d_in += 1
        elif d_in == 2:
            lin_d = line.strip().split(' ')
            if len(lin_d) < 2:
                d_in = 0
                break
            else:
                exclude.append(lin_d[1])
        else:
            continue



## order is exclude, trim, duplicated




fh = pysam.FastxFile(fasta)

out_seqs = set()

for entry in fh:
    if entry.name in c_dict:
        new_seqs = []
        ranges = c_dict[entry.name]
        if len(ranges) > 1:
            for x in ranges:
                start = int(x[0])
                stop = int(x[1])
                new_seq_1 = entry.sequence[:start]  
                new_seqs.append(new_seq_1)
                new_seq_2 = entry.sequence[stop+1:]
                new_seqs.append(new_seq_2)
        else:
            start = int(ranges[0][0])
            stop = int(ranges[0][1])
            new_seq_1 = entry.sequence[:start]  
            new_seqs.append(new_seq_1)
            new_seq_2 = entry.sequence[stop+1:]
            new_seqs.append(new_seq_2)
        count = 0
        for i in new_seqs:
            if len(i) >= 200:
                new_name = re.sub("_length_.*$", '', entry.name) + '_length_' + str(len(i)) + '_trimmed'
                out_seq = '>{0}\n{1}\n'.format(new_name +"_" + str(count), i)
                out_seqs.add(out_seq)
            else:
                continue
            count += 1
            continue
    elif entry.name in exclude:
        continue
    else:
        out_seq = '>{0}\n{1}\n'.format(entry.name, entry.sequence)
        out_seqs.add(out_seq)
        continue
        


    
for i in out_seqs:
    print(i)
