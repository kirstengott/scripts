#!/usr/bin/env python

## check the available cpu percentage and print out how many cores are availbe to use under current usage limitations.
## (number of cores available/2)

import subprocess
import sys

## chech the system average usage every 2 seconds, five times and pull out the average usage over that time frame
t    = subprocess.check_output(["sar", "-u", "2", "5"])
line = t.split("\n")
line.reverse()
ave  = line[1]
new_ave = ave.split()
new_ave.reverse()
perc_cpu_avail = new_ave[0]

## check the total number of cores available on the machine, then calculate what percentage of that you can use.
nproc = subprocess.check_output(["nproc"])
cpu_total = nproc.split('\n')[0]
avail_for_use = str(int((float(perc_cpu_avail)/100 * int(cpu_total))/2))


## write out the number without a new line character (so this can be pulled in by outside programs as a string and be ready to use).
sys.stdout.write(avail_for_use)

