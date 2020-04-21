#!/usr/bin/env python3

import os
import sys
import subprocess

def unzip_file(fasta):
    command = 'gunzip {}'.format(fasta)
    print('Running:', command)
    p = subprocess.Popen(command, shell = True)
    os.waitpid(p.pid, 0)

def run_busco(fasta, threads, run_type, lineage):
    out, ext = os.path.splitext(fasta)
    if ext == '.gz':
        unzip_file(fasta)
        fasta = out
        out, ext = os.path.splitext(out)
        
    if run_type == 'prot':
        command = "docker run -u $(id -u) -v $(pwd):/busco_wd ezlabgva/busco:v4.0.beta1_cv1 run_BUSCO.py -c {} -i {} -o {} -m prot -l {}".format(threads, fasta, out, lineage)
    elif run_type == 'tran':
        command = "docker run -u $(id -u) -v $(pwd):/busco_wd ezlabgva/busco:v4.0.beta1_cv1 run_BUSCO.py -c {} -i {} -o {} -m tran -l {}".format(threads, fasta, out, lineage)
    elif run_type == 'geno':
        command = "docker run -u $(id -u) -v $(pwd):/busco_wd ezlabgva/busco:v4.0.beta1_cv1 run_BUSCO.py -c {} -i {} -o {} -m geno -l {}".format(threads, fasta, out, lineage)
    else:
        print("BUSCO run type not provided, exiting")
        sys.exit(1)
    print('Running:', command)
    p = subprocess.Popen(command, shell = True)
    os.waitpid(p.pid, 0)

def main(argv):
    if len(argv) < 3:
        sys.stderr.write('''Usage: %s <fastafile> <threads> <run_type>
        Runs protein BUSCO from docker container and will unzip fastafile for you
        <fastafile>: input fasta file
        <threads>: number of threads to run on
        <run_type>: BUSCO input option, one of [prot, tran, geno]
        <lineage>: BUSCO lineage to use \n''' % os.path.basename(sys.argv[0]))
        sys.exit(1)
    fasta   = sys.argv[1]
    threads = sys.argv[2]
    if not sys.argv[4]:
        lineage = 'ascomycota_odb10'
    else:
        lineage = sys.argv[4]
    run_busco(fasta = fasta, threads = threads, run_type = sys.argv[3], lineage = lineage)

if __name__ == "__main__":
    main(sys.argv[1:])

