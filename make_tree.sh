#!/bin/bash

usage="Usage: ./make_tree.sh -f <align.faa>
                -h: this help message.
                -f:  folder containing files with amino acid sequences of orthologous genes [REQUIRED]
                -p:  number of threads to use [REQUIRED]
                "


while getopts ":hf:p:" opt; do
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



## setting up directories
mkdir tmp
mkdir astral
mkdir astral/trees
mkdir astral/bootstraps


files=`ls $database`



for file in $files
do

## create basename for file creation
file_strip=`basename $file`
filename=${file_strip%.*}

#echo $file_strip

echo "Running mafft alignment: mafft ${database}/$file >tmp/${filename}.afa "

mafft ${database}/$file >tmp/${filename}.afa

done




## count the number of mafft alignments that are total matches
## grep - -L tmp/* | wc -l

## remove any alignments that are all matches
grep - -L tmp/* | parallel -j 3 rm {}

cd tmp

for file in $files
do

## create basename for file creation
file_strip=`basename $file`
filename=${file_strip%.*}

## convert mafft alignments to phylip format

echo "Converting to maaft alignments to phylip format: aln_convert.pl fasta phylip < tmp/${filename}.afa > tmp/${filename}.phylip "
aln_convert.pl fasta phylip < ${filename}.afa > ${filename}.phylip

## run raxml


echo "Running raxml part 1: raxmlHPC -m GTRGAMMA -n ${filename} -s ${filename}.phylip -f a -x 897543 -N 100 -p 345232 -T ${threads} "

raxmlHPC -m GTRGAMMA -n ${filename} -s ${filename}.phylip -f a -x 897543 -N 100 -p 345232 -T ${threads}


echo "Running raxml part 2: raxmlHPC -f b -m GTRGAMMA -z RAxML_bootstrap.${filename} -t RAxML_bestTree.${filename} -n ${filename}.BS_TREE -T ${threads} "
raxmlHPC -f b -m GTRGAMMA -z RAxML_bootstrap.${filename} -t RAxML_bestTree.${filename} -n ${filename}.BS_TREE -T ${threads}

done

cd ..

echo "Setting up astral file hierarchy "

cd astral/bootstraps
ln -s ../../tmp/RAxML_bootstrap.* ./

cd ../trees
ln -s ../../tmp/RAxML_bipartitions.*BS_TREE ./


## return to the astral root directory
cd ..

echo "Initializing astral run files "

ls bootstraps/* >bs_files
cat trees/* >allgenes.tre



num_genes=`ls trees | wc -l`

echo "Running astral: java -jar ~/bin/astral.5.6.1.jar -i allgenes.tre -b bs_files -o ${num_genes}_genes.tre 2>${num_genes}_genes.log"
java -jar ~/bin/astral.5.6.1.jar -i allgenes.tre -b bs_files -o ${num_genes}_genes.tre 2>${num_genes}_genes.log
tail -n 1 ${num_genes}_genes.tre >bootstrapped_tree_final.tre

echo 'Final consensus tree found at astral/bootstrapped_tree_final.tre'
