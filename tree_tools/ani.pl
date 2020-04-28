#!/bin/env perl

## this script was from Marc Chevrette, I just added in the column with mean ANI
use strict;
use warnings;
use List::Util qw(sum);


my @genomes = (
    @ARGV
    );
open my $ofh, '>', 'ANIvis.tsv' or die $!;
print $ofh join("\t", 'Query', 'Reference', 'ANI', 'QueryAligned')."\n";
foreach my $q (@genomes){
    foreach my $ref (@genomes){
	next if($q eq $ref);
	system("fastANI -q $q -r $ref --visualize -t 10 -o tmp.txt");
	open my $vfh, '<', 'tmp.txt.visual' or die $!;
	my $alen = 0;
	my $ani_array = 0;
	my $ani_length = 0;
	while(<$vfh>){
	    chomp;
	    my @ln = split(/\t/, $_);
	    $alen += abs($ln[6] - $ln[7])+1;
	    $ani_array += $ln[2];
	    $ani_length += 1
	}
	close $vfh;
	my $ani_mean = 0;
	$ani_mean = sum($ani_array)/$ani_length;
	print $ofh join("\t", $q, $ref, $ani_mean, $alen)."\n";
	system("rm tmp.txt*");
    }
}
close $ofh;
