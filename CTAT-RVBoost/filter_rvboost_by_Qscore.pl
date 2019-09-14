#!/usr/bin/env perl

use strict;
use warnings;


my $usage = "\n\tusage: $0 rvboost.vcf [minQscore=0.05]\n\n";

my $rvboost_vcf_file = $ARGV[0] or die $usage;
my $min_qscore = $ARGV[1];

unless (defined $min_qscore) {
    $min_qscore = 0.05;
}


main: {
    
    open(my $fh, $rvboost_vcf_file) or die "Error, cannot open file $rvboost_vcf_file";
    
    while(<$fh>) {
        if (/^\#/) {
            print;
            next;
        }
        /QScore=([\d\.]+)/ or die "Error, no QScore found for $_";
        my $qscore = $1;
        if ($qscore >= $min_qscore) {
            print;
        }
    }

    exit(0);
}


        
        
    
