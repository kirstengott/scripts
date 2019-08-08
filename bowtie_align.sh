#!/bin/bash


usage="Usage: ./bowtie_align.sh -r <ref_fasta> -1 <1> -2 <2> -s
                -h: this help message.
                -r: reference fasta file [REQUIRED]
                -1: forward reads
                -2: reverse reads
                -t: threads
                -s: skip bowtie build flag
                -d: skip samtools depth step
"


while getopts ":hr:1:2:t:sd" opt; do
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

if [ -z "$skip" ]
then
    bowtie2-build $fasta $fasta
else
    echo 'skipping bowtie build index'
fi


if [ -z "$fastq2" ]
then
    bowtie2 -x $fasta -p $threads -U $fastq1 --met-file bt2_alignment_metric.txt | samtools view -b -@ $thread -o $out.bam    
else
    bowtie2 -x $fasta -p $threads -1 $fastq1 -2 $fastq2 --met-file bt2_alignment_metric.txt | samtools view -b -@ $threads -o $out.bam
    fi


samtools sort -n -o ${out}.sorted.bam -T ${out}.tmp ${out}.bam
samtools index ${out}.sorted.bam

if [ -z "$skip" ]
then
samtools depth ${out}.sorted.bam >${out}.depth
gzip ${out}.depth 
else
    echo 'skipping depth estimations'
fi
