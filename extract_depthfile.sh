#!/bin/bash

#Maintainer: Ian Rambo :: ian.rambo@utexas.edu
#Thirteen... that's a mighty unlucky number... for somebody!

usage="$(basename "$0"): extract coverage profiles for a genome.

If extracting for mmgenome2, make sure --out ends with _cov

extract_depthfile.sh -f <FASTA> -d <DEPTHFILE> -v <1|0> > out_cov

where:

   -h --- show this help message
   -f --- FASTA file to extract depths for
   -d --- depth file from jgi_summarize_bam_contig_depths
   -o --- output depthfile subset
   -v --- exclude depthfile variance columns (0|1). Default=1

    "

novar=1

while getopts ':hf:d:o:v:' option; do
    case "${option}" in
    h) echo "$usage"
       exit ;;
    f) fasta=${OPTARG};;
    d) depth=${OPTARG};;
    o) outdepth=${OPTARG};;
    v) novar=${OPTARG};;

    :) printf "missing argument for -%s\n" "$OPTARG" >&2
       echo "$usage" >&2
       exit 1 ;;
    \?) printf "illegal option: -%s\n" "$OPTARG" >&2
        echo "$usage" >&2
        exit 1 ;;
    esac
done

shift $((OPTIND - 1))
#=============================================================================
function check_file() {
    if [ ! -f "$1" ]; then
        echo "could not find $1, no such file"
        exit 1
    else
        return 0
    fi
}
#-----------------------------------------------------------------------------
check_file $fasta && check_file $depth

depthfile_ncol=$(awk -F'\t' '{print NF; exit}' ${depth})

if [ $novar == 1 ]; then
    depth_colseq=$(seq 4 2 ${depthfile_ncol} | tr '\r\n' ',' | sed -e 's/,$//g')
elif [ $novar == 0 ]; then
    depth_colseq="4-"
else
    echo "ERROR: --var must be set to 1 or 0"
    exit 1
fi

#Write the header line (contig ID and sample depths)
head -n 1 $depth | cut -f1,${depth_colseq} -d $'\t' > $outdepth

#Write the depthfile rows for the contigs
grep -f <(grep '>' $fasta | sed -e "s/^>//g" | \
awk '{gsub("$","[[:space:]]*",$0); print;}') $depth | \
cut -f1,${depth_colseq} -d $'\t' >> $outdepth
