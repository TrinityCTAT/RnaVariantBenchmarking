#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;

my $usage = "\n\tusage: $0 gatk_annot_vcf.\n\n";

my $gatk_annot_vcf = $ARGV[0] or die $usage;

main: {

    open(my $fh, $gatk_annot_vcf) or die "Error, cannot open file: $gatk_annot_vcf";
    while(<$fh>) {
        if (/^\#/) { 
            # vcf header
            print; 
            next; 
        }
        chomp;
        my @x = split(/\t/);
        my $annot_field = $x[7];
        my @fields = split(/;/, $annot_field);
        my %keyvals = map { my ($key, $val) = split(/=/, $_); $key => $val } @fields;
        
        #print Dumper(\%keyvals);
        
        ## make adjustments.
        $keyvals{RPT} = (defined $keyvals{RPT}) ? 1:0;
        $keyvals{SPLICEADJ} = (defined $keyvals{SPLICEADJ}) ? $keyvals{SPLICEADJ} : -1;
        
        my @new_vals;
        foreach my $key (keys %keyvals) {
            my $val = $keyvals{$key};
            unless (defined $val) {
                #print STDERR "Warning, missing value for $key in: $annot_field\n";
                next;
            }
            push (@new_vals, "$key=$val");
        }
        $x[7] = join(";", sort @new_vals);

        print join("\t", @x) . "\n";
                
    }

    exit(0);
}


