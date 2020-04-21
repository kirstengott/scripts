#!/bin/bash


usage="Usage: ./get_sra.sh -f <SRR file> -d <outdir>
                -h: this help message.
                -f: File with SRR accessions [REQUIRED]
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

mkdir $dir

for id in $ids
do
    base=`echo $id | sed -E "s/(SRR[0-9]{3}).*$/\1/"`
    wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/${base}/${id}/${id}.sra
    fastq-dump --outdir $dir --gzip --skip-technical --readids --read-filter pass --dumpbase --split-3 --clip ${id}.sra #dump the fastqs for the fetched SRA
    #prefetch $id #prefetch the SRA id alternative (slow)
    #fastq-dump --outdir $dir --gzip --skip-technical --readids --read-filter pass --dumpbase --split-3 --clip ${id} #dump the fastqs for the fetched SRA
done
