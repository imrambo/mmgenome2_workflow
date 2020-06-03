#!/bin/bash

#Maintainer: Ian Rambo :: ian.rambo@utexas.edu
#Thirteen... that's a mighty unlucky number... for somebody!


usage="$(basename "$0"): extract coverage profiles for a genome.

If extracting for mmgenome2, make sure --out ends with _cov

extract_depthfile.sh -f <FASTA> -d <DEPTHFILE> -v <1|0> > out_cov

where:

   -h | --help  --- show this help message
   -f | --fasta --- FASTA file to extract depths for
   -d | --depth --- depth file from jgi_summarize_bam_contig_depths
   -o | --out   --- output depthfile subset
   -v | --novar --- exclude depthfile variance columns (0|1). Default=1

    "

novar=1

while [ "$1" != "" ]; do
    case $1 in
        -f | --fasta )          shift
                                fasta=$1
                                ;;
        -d | --depth )          shift
                                depth=$1
                                ;;
        -o | --out )            shift
                                outdepth=$1
                                ;;
        -v | --var )            shift
                                novar=$1
                                ;;
        -h | --help )           echo "$usage"
                                exit
                                ;;
        * )                     echo $usage
                                exit 1
    esac
    shift
done
#=============================================================================

depthfile_ncol=$(awk -F'\t' '{print NF; exit}' ${depth})

if [ $novar == 1 ]; then
    depth_colseq=$(seq 4 2 ${depthfile_ncol} | tr '\r\n' ',' | sed -e 's/,$//g')
elif [ $novar == 0 ]; then
    depth_colseq="4-"
else
    echo "ERROR: --var must be set to 1 or 0"
    exit 1
fi

#Write the output depthfile header
head -n 1 $depth | cut -f1,${depth_colseq} -d $'\t' > $outdepth
#Write the depthfile columns for the contigs
grep '>' $fasta | sed -e "s/^>//g" | \
awk '{gsub("$","[[:space:]]*",$0); print;}' | \
grep -f - $depth | \
cut -f1,${depth_colseq} -d $'\t' >> $outdepth
