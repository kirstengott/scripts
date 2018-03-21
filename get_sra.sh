#!/bin/bash


usage="Usage: ./get_sra.sh -f <SraRunInfo.csv>
                -h: this help message.
                -f: SraRunInfo.csv [REQUIRED]"


while getopts ":hf:" opt; do
    case $opt in
	h)
	    echo "$usage" >&2
	    exit 1
	    ;;
	f)
	    file="$OPTARG" >&2
	    echo "SraRunInfo.csv file: $OPTARG" >&2
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

ids=`cut -d , -f 1 $file | grep -v Run`

for id in $ids
do
    prefetch $id #prefetch the SRA id
    fastq-dump --split-3 --gzip $id #dump the fastqs for the fetched SRA
done
