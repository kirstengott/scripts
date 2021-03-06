import sys, subprocess
def run_trimal(afa):
    ''' Runs TrimAl to trim the alignments. 
    afa: the output from aligning with mafft
    recompute: whether or not to recopute the trimming
    '''
    output_file = afa + '.trimmed'
    if os.path.exists(output_file):
        pass
    else:
        command = "trimal -in {} -out {} -automated1".format(afa, output_file)
        print('Running:', command)
        p = subprocess.Popen(command, shell=True)
        pid = os.waitpid(p.pid, 0)
        pid = list(pid)[1]
        if pid != 0:
            print('trimal failed, check input fasta file')
            sys.exit()
    return(output_file)


run_trimal(sys.argv[1])

