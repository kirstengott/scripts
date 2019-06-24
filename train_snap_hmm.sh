#!/bin/bash

usage="Usage: ./train_snap_hmm.sh -p <prefix> -d <*_master_datastore_index.log>
                -h: this help message.
                -p: output prefix
                -d: maker *_master_datastore_index.log [REQUIRED]
"


while getopts ":hp:d:" opt; do
    case $opt in
	h)
	    echo "$usage" >&2
	    exit 1
	    ;;
	p)
	    pref="$OPTARG" >&2
	    ;;
	d)
	    db="$OPTARG" >&2
	    ;;
	\?)
	    echo "Invalid option: -$OPTARG" >&2
	    echo "$usage" >&2
	    exit 1
	    ;;
	:)
	    echo "Option -$OPTARG requires an argument." >&2
	    echo "$usage" >&2
	    exit 1
	    ;;
    esac
done



maker2zff -c 0 -e 0 -l 50 -d $db
fathom genome.ann genome.dna -gene-stats >gene-stats.log
fathom genome.ann genome.dna -validate >validate.log 2>&1
fathom genome.ann genome.dna -categorize 1000 > categorize.log 2>&1
fathom uni.ann uni.dna -export 1000 -plus > uni-plus.log 2>&1
mkdir params; cd params
forge ../export.ann ../export.dna > ../forge.log 2>&1
hmm-assembler.pl $pref params >${pref}.hmm
