#!/bin/bash

usage="Usage: ./split_fasta.sh -f <fastafile>
                Splits a FASTA file into individual files for each record.
                -h: this help message.
                -f: input fasta to split [REQUIRED]
"


while getopts ":hf:" opt; do
    case $opt in
	h)
	    echo "$usage" >&2
	    exit 1
	    ;;
	f)
	    myseq="$OPTARG" >&2
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


while read line
do
    if [[ ${line:0:1} == '>' ]]
    then
	outfile=${line#>}.fa
	echo $line > $outfile
    else
	echo $line >> $outfile
    fi
    done < $myseq
