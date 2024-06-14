#!/usr/bin/bash
# Author: Jinxin Meng
# Created Date: 2023-09-30

if [ $# -ne 3 ];then
    echo "$0 <fq1,fq2> <index> <out_prefix>"
    exit 127
fi

fq=$1
fq1=$(echo $fq | cut -d "," -f1)
fq2=$(echo $fq | cut -d "," -f2)
idx=$2
p=$3
trds=16

bowtie2 --end-to-end --mm --fast -p $trds -1 $fq1 -2 $fq2 -x $idx 2> $p.log |\
    samtools view -@ $trds -bS | samtools sort -@ $trds -o ${p}_sort.bam - 2>/dev/null
chmod 444 $p.log ${p}_sort.bam
