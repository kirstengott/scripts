#!/bin/bash


usage="Usage: ./run_targetp.sh -f <fastafile>
                Splits a FASTA file into individual files for each record.
                -h: this help message.
                -f: input fasta [REQUIRED]
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


rm -rf tp_tmp
mkdir tp_tmp
mkdir output

base=`basename $myseq`
echo $base

cd tp_tmp
rename_fasta_seqs.py ../$myseq $base ${base}_idmap.txt
split_fasta_to_sizes.py ${base} 1000
rm $base
cd ..



for file in `ls tp_tmp/sequences*fa`
do
    file_h=`basename $file`
    ./bin/targetp -N tp_tmp/$file_h >tp_tmp/${base}_$file_h
    rm tp_tmp/$file_h
done


echo -e 'Name\tLen\tmTP\tSP\tother\tLoc\tRC' >>output/$base

for file in `ls tp_tmp/*fa`
do
    file_h=`basename $file`
    ./parse_targetp.py tp_tmp/$file_h tp_tmp/${base}_idmap.txt >>output/$base
done

rm -rf tp_tmp
