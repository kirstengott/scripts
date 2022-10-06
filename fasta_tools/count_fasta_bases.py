import pysam, sys, os


def get_fasta_length(fasta):
    fh = pysam.FastxFile(fasta)
    total_length = 0
    for entry in fh:
        total_length += len(entry.sequence)
    return(total_length)


def main(argv):
    if len(argv) != 1:
        sys.stderr.write("Usage: %s <fastafile> \n" % os.path.basename(sys.argv[0]))
        sys.exit(1)
    length = get_fasta_length(fasta = sys.argv[1])
    print(sys.argv[1], length)


if __name__ == '__main__':
    main(sys.argv[1:])    
        
