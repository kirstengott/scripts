#!/bin/bash

## paths to executables

mafft="/usr/bin/mafft"
raxml="/usr/bin/raxmlHPC"
catfasta="/home/kgotting/scripts/catfasta2phyml.pl"
aln_convert="/home/kgotting/scripts/aln_convert.pl"

usage="Usage: ./make_tree.sh -f <align.faa>
                -h: this help message.
                -f:  folder containing sequence files to align [REQUIRED]
                -c:  number of threads to use [REQUIRED]
                -a: skip mafft alignment step
                -p: input is protein sequences
                "


while getopts ":hf:pc:a" opt; do
    case $opt in
	h)
	    echo "$usage" >&2
	    exit 1
	    ;;
	f)
	    database="$OPTARG" >&2
	    echo "Database: $OPTARG" >&2
	    ;;
	c)
	    threads="$OPTARG" >&2
	    echo "Number of threads: $OPTARG" >&2
	    ;;
	a)
	    align=1
	    ;;
	p)
	    protein=1
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

mkdir tmp

if [ -z "$align" ]
then
    for file in $files
    do
	## create basename for file creation
	file_strip=`basename $file`
	filename=${file_strip%.*}
	echo "Running mafft alignment: ${mafft} ${database}/$file >tmp/${filename}.afa "
	${mafft} ${database}/$file >tmp/${filename}.afa

    done

    ## remove any alignments that are all matches
    grep - -L tmp/* | parallel -j 3 rm {}


else
    echo 'skipping maaft alignment'
fi



cd tmp

## cat into one fasta and convert it into phylip format
${catfasta} -f *afa >../all.afa
cd ..

${aln_convert} fasta phylip < all.afa > all.phylip



if [ -z "$protein" ]; then
    echo 'Running RAXML on nucleotides'
    ## I want to keep all of the leaves, even if everything matches.
    ${raxml} -m GTRGAMMA -n all -s all.phylip -f a -x 897543 -N autoMRE -p 345232 -T ${threads}
    ${raxml} -f b -m GTRGAMMA -z RAxML_bootstrap.all -t  RAxML_bestTree.all -n all.BS_TREE -T ${threads}
else
    ## proteins
    echo 'Running RAXML on proteins'
    ## -f a rapid Bootstrap analysis and search for best­scoring ML tree in one program run
    ${raxml} -m PROTGAMMAAUTO -n all -s all.phylip -f a -x 897543 -N autoMRE -p 345232 -T ${threads}
    ${raxml} -f b -m PROTGAMMAAUTO -z RAxML_bootstrap.all -t  RAxML_bestTree.all -n all.BS_TREE -T ${threads}

fi
