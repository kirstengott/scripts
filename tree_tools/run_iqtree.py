import sys, subprocess

def run_iqtree(phylip, threads):
    ''' Runs iqtree with phylip file, and number of threads provided.
    '''
    if not os.path.exists('iqtrees'):
        os.mkdir('iqtrees')
    output_dir = os.path.join(os.getcwd(), 'iqtrees')
    name = re.sub('.phylip', '', os.path.basename(phylip))
    prefix = output_dir + "/" + name
    output_filename = name + '.treefile'
    output_filename = os.path.join(output_dir, output_filename)
    if os.path.exists(output_filename):
        pass
    else:
        command1 = "iqtree -s {phylip} -nt {threads} -bb 1000 -alrt 1000 -pre {outp} > {outp}_stdout_iqtree.txt".format(phylip = phylip, threads = threads, outp = prefix)
        print('iqtree command:', command1)
        p = subprocess.Popen(command1, shell = True)
        os.waitpid(p.pid, 0)
    return(output_filename)
