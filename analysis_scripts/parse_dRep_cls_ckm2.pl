#!/usr/bin/perl
# Author: Jinxin Meng
# Create Date: 2023-04-04
# Modified Date: 2023-11-07
# parse dRep cluster result, and a representative genome is selected if it has the 
# highest quality score (completeness - 5 * contanmination) provided by checkm2.
# dRep perform two clsuter, the result file is different from previous file. Thus,
# in this version, we amend the some code.


use warnings;
use strict;

die "Usage: perl $0 [Cdb.csv] [checkm2 quality_report.tsv] [out_file]\n" unless @ARGV eq 3;
my ($cdb, $ckm, $out) = @ARGV;

open I,"<$cdb" or die "Can't open $cdb: $!\n";
open IN,"<$ckm" or die "Can't open $ckm: $!\n";
open O, ">$out";

my (%cls, %qs, @genome) = (); 
while (<I>) {
    chomp;
    next if /primary_cluster/;
    my @s = split /,/;
    $s[0] =~ s/.(fa|fna)//;
    push @{$cls{$s[1]}}, $s[0];
    push @genome, $s[0];
}

while (<IN>) {
    chomp;
    next if /Completeness/;
    my @s = split /\t/;
    $qs{$s[0]} = $s[1] - 5 * $s[2] if (grep {$_ eq $s[0]} @genome);
}

for my $i (sort keys %cls) {
    if (@{$cls{$i}} eq 1) {
        print O "cluster_$i\t1\t@{$cls{$i}}\t@{$cls{$i}}\n";
    } else {
        my $len = 0;
        my $max = 0;
        my $g = ();
        for my $j (@{$cls{$i}}) {
            if ($qs{$j} >= $max) {
                $max = $qs{$j};
                $g = $j;
                $len+=1;
            } else {
                $len+=1;
            }
        }
        print O "cluster_$i\t$len\t$g\t".join(",", @{$cls{$i}})."\n";     
    }   
}
close(O);
close(I);
close(IN);

