
import sys

line1 = 0
with open(sys.argv[1], 'r') as fh:
    for line in fh:
        if line1 == 0:
            line1 += 1
        else:
            line_in = line.strip().split()
            print(">{}\n{}\n".format(line_in[0], line_in[1]))
        
