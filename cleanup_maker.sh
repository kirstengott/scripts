#!/bin/bash

usage="Usage: ./cleanup_maker.sh -d <maker output directory> -o <path to output directory>
                -h: this help message.
                -d: the maker-created output directory [REQUIRED]
                -o: the desired output location"


while getopts ":hd:o:" opt; do
    case $opt in
	h)
	    echo "$usage" >&2
	    exit 1
	    ;;
	d)
	    i="$OPTARG" >&2
	    echo "Maker directory: $OPTARG" >&2
	    ;;
	o)
	    o="$OPTARG" >&2
	    echo "Output file path: $OPTARG" >&2
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


cd $i

gff3_merge -d *master_datastore_index.log -g -n  ## pull out the gene models
evi=`ls *all.gff`
mv $evi ${evi%all.gff}genemodels.gff
gff3_merge -d *master_datastore_index.log -n ## pull out the evidence
fasta_merge -d *master_datastore_index.log; ## pull out the fasta files


if [ ! -d "$o" ]; then
    mkdir $o
fi



## copy the files to an output directory
cp *all.maker.proteins.fasta $o
cp *all.maker.transcripts.fasta $o
cp *all.gff $o
cp *.genemodels.gff $o

