#!/usr/bin/perl
# encoding: utf-8
# author: Jinxin Meng
# created date: 2023-12-21, 19:42:50
# modified date: 2023-12-21, 20:45:02
use warnings;
use strict;

die "Usage: perl $0 [gff] [out_file]\n" unless @ARGV eq 2;

my ($in_f, $out_f) = @ARGV;
open I, "<$in_f" or die "Cannot open $in_f: $!\n";
open O, ">$out_f";
my ($flag, $x) = "";
while (<I>) {
    next if /^#/;
    my @s = split/\t/;
    if ($s[0] ne $flag) {
        $x = 1;
        $flag = $s[0];
    } else {
        $x++;
    }
    print O "$s[0]\t$s[0]_$x\t$s[3]\t$s[4]\t$s[6]\n"; 
}

