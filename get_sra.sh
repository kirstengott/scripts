#!/bin/bash


usage="Usage: ./get_sra.sh -f <SraRunInfo.csv>
                -h: this help message.
                -f: SraRunInfo.csv [REQUIRED]
                -d: outdir [required]

"


while getopts ":hf:d:" opt; do
    case $opt in
	h)
	    echo "$usage" >&2
	    exit 1
	    ;;
	f)
	    file="$OPTARG" >&2
	    echo "SRR accession file: $OPTARG" >&2
	    ;;
	 d)
	    dir="$OPTARG" >&2
	    echo "SRR accession file: $OPTARG" >&2
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

ids=`cat $file`

for id in $ids
do
    first_ind=`echo $id | grep -o -P 'SRR[0-9]{3}'`
    wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/${first_ind}/${id}/${id}.sra
    #prefetch $id #prefetch the SRA id
    fastq-dump --outdir $dir --gzip --skip-technical --readids --read-filter pass --dumpbase --split-3 --clip $id #dump the fastqs for the fetched SRA
done
