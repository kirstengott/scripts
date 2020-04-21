#!/bin/bash


usage="Usage: ./bowtie_align.sh -i <ref_id_file> -o <output_dir>
                -h: this help message.
                -i: unique ids for files to process
                -o: output directory
                -p: number of threads
                -g: path to genomescope bin
"


while getopts ":hi:o:p:g:" opt; do
    case $opt in
	h)
	    echo "$usage" >&2
	    exit 1
	    ;;
	i)
	    file="$OPTARG" >&2
	    ;;
	o)
	    out="$OPTARG" >&2
	    ;;
	p)
	    threads="$OPTARG" >&2
	    ;;
	g)
	    gs="$OPTARG" >&2
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


mkdir $out

for i in `cat $file`;
do echo $i
   ls $i*fastq.gz | parallel gunzip {}
   jellyfish count -C -m 31 -s 3G -t $threads *.fastq -o ${out}/${i}.reads.jf;
   kat hist ${out}/${i}.reads.jf -t $threads -o ${out}/$i.histo
   kat gcp -t $threads -o ${out}/$i.histo ${out}/${i}.reads.jf
   ls $i*fastq | parallel gzip {}
   grep -v "#" $i >${i}.fix
   $gs -i $i -o ${i%histo.fix}genomescope -k 31
   rm ${out}/${i}.reads.jf
done
