#!/usr/bin/env python3

import os
import sys
import subprocess


def run_busco(fasta, threads):
    out, ext = os.path.splitext(fasta)
    out = out + "_busco"
    command = "docker run -u $(id -u) -v $(pwd):/busco_wd ezlabgva/busco:v4.0.beta1_cv1 run_BUSCO.py -c {} -i {} -o {} -m prot -l ascomycota_odb10".format(threads, fasta, out)
    print('Running:', command)
    p = subprocess.Popen(command, shell = True)
    os.waitpid(p.pid, 0)

def main(argv):
    if len(argv) != 2:
        sys.stderr.write("Usage: %s <fastafile> <threads>\n" % os.path.basename(sys.argv[0]))
        sys.exit(1)
    fasta   = sys.argv[1]
    threads = sys.argv[2]
    run_busco(fasta = fasta, threads = threads)

if __name__ == "__main__":
    main(sys.argv[1:])

