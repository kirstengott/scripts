#!/bin/bash

usage="Usage: ./spades_pipe.sh -p <prefix> -1 <1> -2 <2> -u <u>
                -h: this help message.
                -p: prefix for output file [REQUIRED]
                -1: forward fastq
                -2: reverse fastq
                -u: unpaired reads
"


while getopts ":hp:1:2:u:" opt; do
    case $opt in
	h)
	    echo "$usage" >&2
	    exit 1
	    ;;
	p)
	    pref="$OPTARG" >&2
	    ;;
	1)
	    fwd="$OPTARG" >&2
	    ;;
	2)
	    rev="$OPTARG" >&2
	    ;;
	u)
	    unm="$OPTARG" >&2
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


k=31
threads=10


gunzip $fwd
gunzip $rev


echo 'Running musket\n'
musket -p $threads -inorder -omulti $pref.musket -k $k 536870912 $fwd $rev

mv $pref.musket.0 $pref.musket.0.fastq
mv $pref.musket.1 $pref.musket.1.fastq


flash -t $threads -o $pref $pref.musket.0.fastq $pref.musket.1.fastq

gunzip $unm
cat $pref.extendedFrags.fastq ${unm%.gz} >$pref.unmated_reads.fastq
gzip $unm

spades.py --pe1-1 $pref.notCombined_1.fastq --pe1-2 $pref.notCombined_2.fastq --s2 $pref.unmated_reads.fastq -t $threads -o $pref-spades

gzip *fastq

cp $pref-spades/contigs.fasta ./$pref-spades.fna

## align the reads to the new assembly


## run bowtie2
echo 'Aligning DNA reads to genome assembly'


bowtie_align.sh -r ./$pref-spades.fna -1 $fwd -2 $rev -t 10


## build the bowtie index
#bowtie2-build $pref-spades.fna $pref-spades.fna

# bowtie2 -x ./$pref-spades.fna -p 10 -1 $fwd -2 $rev | samtools view -b -@ 10 -o $pref.bam

# samtools sort -o $pref.sorted.bam $pref.bam
# samtools index $pref.sorted.bam

# echo 'Calculating Depth of Coverage'
# samtools depth $pref.sorted.bam >$pref.depth

# gzip $pref.depth







