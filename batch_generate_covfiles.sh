#!/bin/bash 

# Maintainer: Ian Rambo :: ian.rambo@utexas.edu
# Thirteen... that's a mighty unlucky number... for somebody!

usage="$(basename "$0"): Parallel batch script to generate coverage files for genomic bins meeting certain quality threshold criteria

where:

   -h --- show this help message
   -c --- CheckM output file with Completeness and Contamination
   -b --- input FASTA bin directory 
   -d --- depth file from jgi_summarize_bam_contig_depths
   -e --- FASTA extension. Default=fa
   -o --- output base directory
   -v --- exclude depthfile variance columns (0|1). Default=1
   -m --- genome completeness minimum
   -n --- genome completeness maximum
   -y --- genome contamination minimum
   -z --- genome contamination maximum
   -j --- number of parallel jobs to execute. Default=2
   -l --- directory to write GNU parallel joblog to

    "
# Batch script to generate coverage files for genomic bins 
# meeting certain quality threshold criteria.

has_command () {
    command -v "$1" >/dev/null 2>&1 || \
    { echo "Requires $1. Ensure that $1 is in your \$PATH."; exit 1; }
}

has_command parallel

pjobs=2
ext="fa"
novar=1
joblog_dir=$(pwd)

while getopts ':hf:d:o:v:' option; do
    case "${option}" in
    h) echo "$usage"
       exit ;;
    c) checkm_file=${OPTARG};;
    b) bindir=${OPTARG};;
    d) depthfile=${OPTARG};;
    e) ext=${OPTARG};;
    v) novar=${OPTARG};;
    o) outdir=${OPTARG};;
    j) pjobs=${OPTARG};;
    l) joblog_dir=${OPTARG};;
    m) comp_min=${OPTARG};;
    n) comp_max=${OPTARG};;
    y) cont_min=${OPTARG};;
    z) comp_max=${OPTARG};;
    :) printf "missing argument for -%s\n" "$OPTARG" >&2
       echo "$usage" >&2
       exit 1 ;;
    \?) printf "illegal option: -%s\n" "$OPTARG" >&2
        echo "$usage" >&2
        exit 1 ;;
    esac
done

shift $((OPTIND - 1))

test -f $checkm_file || echo "CheckM file $checkm_file not found"
test -d $bindir || echo "Directory $bindir containing input genomic bins not found"
test -f $depthfile || echo "Input depth file $depthfile not found"

tstamp=$(date +'%Y-%m-%d_%H-%M-%S')

joblog=${joblog_dir}/$(basename "$0")_${tstamp}.joblog


tail -n +4 $checkm_file  | \
    awk -v comp_min="$comp_min" -v comp_max="$comp_max" -v cont_min="$cont_min" -v cont_max="$cont_max" -F '[[:space:]]+' '($14 >= comp_min) && ($14 <= comp_max) && ($15 >= cont_min) && ($15 <= cont_max) {print $2}' | \
    parallel --joblog $joblog --jobs $pjobs test -d ${outdir}/{} '||' mkdir -p ${outdir}/{}';' bash $PWD/extract_covfile.sh -f ${bindir}/{}.${ext} -d $depthfile -o ${outdir}/{}/{}_cov -v 1 
