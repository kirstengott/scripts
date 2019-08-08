#!/usr/bin/env python
import re
import sys
import statistics

## take in output from samtools depth and output average depth across windows of size sys.argv[2]

if sys.argv[1] == '-':
    in_f = sys.stdin
else:
    in_f = open(sys.argv[1])

window     = int(sys.argv[2])

    

print('contig\twindow_length\twindow\tdepth')

contig = ''
depths = []
count = 0
window_count = 0



for line in in_f:
    items = line.strip().split("\t")
    ## make switch contig logic, reset depths and contig val for each new contig
    contig_item = items[0]
    depth = float(items[2])
    if contig_item == contig: ## we have seen this contig and need to increment our count and add depth to list
        if count <= window:
            depths.append(depth)
        if count == window:
            print("{}\t{}\t{}\t{}".format(contig, count, window_count, statistics.mean(depths)))
            window_count += 1 
            depths = []
            count = 1
        count += 1
    else:
        contig = contig_item ## else reset everything and print out the depths list
        count = 0
        if len(depths) != 0:
            print("{}\t{}\t{}\t{}".format(contig, count, window_count, statistics.mean(depths)))
        depths = [depth]
        window_count = 0
        
    
## recover the last line
if len(depths) != 0:
    print("{}\t{}\t{}\t{}".format(contig, count, window_count, statistics.mean(depths)))
    
