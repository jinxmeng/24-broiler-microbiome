#!/usr/bin/bash

if [ $# -lt 3 ];then
    echo "$0 <fq1> <index> <out_prefix>"
    exit 127
fi

fq1=$1
idx=$2
p=$3
bowtie2 --end-to-end --mm --no-head --no-unal --no-sq \
    -u 20000000 -U $fq1 -x $idx -p 16 2> $p.log |\
    perl -e 'while(<>){if(/^\S+\s+\S+\s+(\S+)\s+/){$h{$1}++;}} for(sort keys %h){print "$_\t$h{$_}\n";}' > $p.rc
chmod 444 $p.log $p.rc 
