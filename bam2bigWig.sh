#!/usr/bin/bash


## define required program locations
bedGraphToBigWig='/n/local/bin/bedGraphToBigWig'
bedItemOverlapCount='/n/local/bin/bedItemOverlapCount'
sort='/usr/bin/sort'
awk='/usr/bin/awk'


usage="Usage: ./bam2bigWig -b <input_bamfile> -d <genome_database> -c <genome_chromosome_file> -s
                -h: this help message.
                -b: the bam file to convert. [REQUIRED]
                -d: the genome database. [REQUIRED]
                -c: location of chromosome lengths file. [REQUIRED]
                -s: flag to make bigWigs for the plus and minus strands"


while getopts ":hb:d:c:s" opt; do
    case $opt in
	h)
	    echo "$usage" >&2
	    exit 1
	    ;;
	b)
	    input_bam="$OPTARG" >&2
	    echo "Bam File: $OPTARG" >&2
	    ;;
	d)
	    genome_name="$OPTARG" >&2
	    echo "Database: $OPTARG" >&2
	    ;;
	c)
	    chromosome_file="$OPTARG" >&2
	    echo "Chromosome File: $OPTARG" >&2
	    ;;
        s)
            strand='TRUE'
            echo "Creating biwWigs for both strands" >&2
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


bam_base=${input_bam%.bam}
output_bed=${bam_base}.bed





## I first convert bam files to bed files. 

echo 'Making BED file'

if [[ -e ${output_bed} ]]; then
    echo 'BED file exists, moving to make bedgraph'
else
    bamToBed -i ${input_bam} -split > ${output_bed}
fi



if [[ $strand == "TRUE" ]]; then

    echo 'Making plus bedgraph'
    ${awk} '{if($6=="+") print}' ${output_bed} | ${sort} -k1,1 | ${bedItemOverlapCount} ${genome_name} -chromSize=${chromosome_file} stdin | ${sort} -k1,1 -k2,2n >${bam_base}.plus.bedGraph


    echo "Making plus bigwig"
    ${bedGraphToBigWig} ${bam_base}.plus.bedGraph ${chromosome_file} ${bam_base}.plus.bw


    echo 'Making minus bedgraph'
    ${awk} '{if($6=="-") print}' ${output_bed} | ${sort} -k1,1 | ${bedItemOverlapCount} ${genome_name} -chromSize=${chromosome_file} stdin | ${sort} -k1,1 -k2,2n | ${awk} '{OFS="\t"; print $1,$2,$3,"-"$4}' >${bam_base}.minus.bedGraph

    echo 'Making minus bigwig'
    ${bedGraphToBigWig} ${bam_base}.minus.bedGraph ${chromosome_file} ${bam_base}.minus.bw

else
    echo 'Making bedgraph'

    cat ${output_bed} | ${sort} -k1,1 | ${bedItemOverlapCount} ${genome_name} -chromSize=${chromosome_file} stdin | ${sort} -k1,1 -k2,2n | ${awk} '{OFS="\t"; print $1,$2,$3,"-"$4}' >${bam_base}.bedGraph

    echo 'Making bigwig'

    ${bedGraphToBigWig} ${bam_base}.bedGraph ${chromosome_file} ${bam_base}.bw

fi
    

#cat ${output_bed} | ${sort} -k1,1 | ${bedItemOverlapCount} ${genome_name} -chromSize=${chromosome_file} stdin | ${sort} -k1,1 -k2,2n >${bam_base}.bedGraph


#${bedGraphToBigWig} ${bam_base}.bedGraph ${chromosome_file} ${bam_base}.bw






