#!/bin/bash


usage="Usage: ./bwa_align.sh -r <ref_fasta> -1 <1> -2 <2> -s -d -o <out>
                -h: this help message.
                -r: reference fasta file [REQUIRED]
                -1: forward reads [REQUIRED]
                -2: reverse reads [REQUIRED]
                -t: threads [REQUIRED]
                -o: outfile prefix [REQUIRED]
                -s: skip bwa build index
                -d: skip samtools depth step
"


while getopts ":hr:1:2:t:sdo:" opt; do
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
	t)
	    threads="$OPTARG" >&2
	    ;;
	s)
	    skip=1
	    ;;
	d)
	    depth=1
	    ;;
	o)
	    out="$OPTARG" >&2
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

#out=`basename $fasta | sed -e "s/.fa//"`

if [ -z "$skip" ]
then
    bwa index $fasta
else
    echo 'skipping bwa build index'
fi


if [ -z "$fastq2" ]
then
    bwa mem $fasta $fastq1 -t $threads | samtools view -b -@ $threads -o $out.bam
    #samtools sort -o ${out}.sorted.bam -T ${out}.tmp ${out}.bam
else
    bwa mem $fasta $fastq1 $fastq2 -t $threads | samtools view -b -@ $threads -o $out.bam
    #samtools sort -o ${out}.sorted.bam -T ${out}.tmp ${out}.bam
    fi



if [ -z "$depth" ]
then
    samtools sort -o ${out}.sorted.bam -T ${out}.tmp ${out}.bam
    samtools index ${out}.sorted.bam
    samtools depth ${out}.sorted.bam >${out}.depth
    gzip ${out}.depth 
else
    echo 'skipping depth estimations'
fi
