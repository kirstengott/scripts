#!/bin/bash

## paths to executables

mafft="/usr/bin/mafft"
raxml="/usr/bin/raxmlHPC"
catfasta="/home/kgotting/scripts/catfasta2phyml.pl"
aln_convert="/home/kgotting/scripts/aln_convert.pl"

usage="Usage: ./make_tree.sh -f <align.faa>
                -h: this help message.
                -f:  folder containing sequence files to align [REQUIRED]
                -p:  number of threads to use [REQUIRED]
                "


while getopts ":hf:p:a" opt; do
    case $opt in
	h)
	    echo "$usage" >&2
	    exit 1
	    ;;
	f)
	    database="$OPTARG" >&2
	    echo "Database: $OPTARG" >&2
	    ;;
	p)
	    threads="$OPTARG" >&2
	    echo "Number of threads: $OPTARG" >&2
	    ;;
	a)
	    align=1
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

## exit if directory isn't set
if [ -z "$database" ]; then
    echo 'Required argument not supplied, see ./make_tree.sh -h'
    exit 1
fi
if [ -z "$threads" ]; then
    echo 'Required argument not supplied, see ./make_tree.sh -h'
    exit 1
fi


files=`ls $database`

if [ -z "$align"]
then
    for file in $files
    do
	mkdir tmp
	## create basename for file creation
	file_strip=`basename $file`
	filename=${file_strip%.*}
	echo "Running mafft alignment: ${mafft} ${database}/$file >tmp/${filename}.afa "
	${mafft} ${database}/$file >tmp/${filename}.afa

    done

    ## remove any alignments that are all matches
    grep - -L tmp/* | parallel -j 3 rm {}

    cd tmp

    ## cat into one fasta and convert it into phylip format
    ${catfasta} -f *afa >all.afa
    ${aln_convert} fasta phylip < all.afa > all.phylip

    cd ..

    else:
    echo 'skipping maaft alignment'
fi

## I want to keep all of the leaves, even if everything matches.
${raxml} -m GTRGAMMA -n all -s tmp/all.phylip -f a -x 897543 -N autoMRE -p 345232 -T ${threads}
${raxml} -f b -m GTRGAMMA -z RAxML_bootstrap.all -t  RAxML_bestTree.all -n all.BS_TREE -T ${threads}
