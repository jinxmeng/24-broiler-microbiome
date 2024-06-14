#!/usr/bin/bash
# Author: Jinxin Meng
# Created Date: 2023-09-30

if [ $# -ne 2 ];then
    echo "$0 <sort.bam> <out_prefix>"
    exit 127
fi

bam=$1
p=$2

/share/data1/software/coverm-0.6.1/coverm genome \
    --bam-files $bam \
    --genome-fasta-directory /share/data1/mjx/proj/01.broiler_metagenome_20230301/05.dRep/species_fna/ \
    --methods length count covered_fraction covered_bases tpm \
    --min-covered-fraction 0 \
    --output-file ${p}_cvm
