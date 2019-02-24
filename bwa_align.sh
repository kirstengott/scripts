#!/bin/bash


usage="Usage: ./bwa_align.sh -r <ref_fasta> -1 <1> -2 <2>
                -h: this help message.
                -r: reference fasta [REQUIRED]
                -1: single end fastq.gz or first mate pair [REQUIRED]
                -2: second mate pair [optional]
"


while getopts ":hr:1:2:" opt; do
    case $opt in
	h)
	    echo "$usage" >&2
	    exit 1
	    ;;
	r)
	    fasta="$OPTARG" >&2
	    ;;
	1)
	    fastq1="$OPTARG" >&2
	    ;;
	2)
	    fastq2="$OPTARG" >&2
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

out=`basename $fasta | sed -e "s/.fa//"`


if [ -z "$fastq2" ]
then
    bwa mem -t 1 $fasta $fastq1 | samtools view -b -@ 10 -o ${out}.bam -
else
    bwa mem -t 1 $fasta $fastq1 $fastq2 | samtools view -b -@ 10 -o ${out}.bam -
    fi


samtools sort -o ${out}.sorted.bam -T ${out}.tmp ${out}.bam
rm ${out}.bam
samtools index ${out}.sorted.bam
samtools depth ${out}.sorted.bam >${out}.depth
gzip ${out}.depth 
