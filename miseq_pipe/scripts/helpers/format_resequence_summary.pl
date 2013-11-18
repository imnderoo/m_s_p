#!/usr/bin/perl -w

use strict;
use warnings;

if ($#ARGV != 1)
{
	print "Usage format_reseqeunce_summary.pl <coverage.summary> <depth of coverage threshold>\n";
	print "$#ARGV parameters found";
	exit;
}

my $coverage = $ARGV[0];
my $threshold = $ARGV[1];

open (COVERAGEFILE, $coverage) or die($!);

my $failedexons = "";
my $passedexons = "";

 
while(my $line = <COVERAGEFILE>) {

	chomp($line);
	
	my @col = split('\t', $line);

	my $min=$col[10];

	if($min ne "min")
	{
		if($min < $threshold)
		{
			$failedexons = $failedexons . "FAIL,$col[0],$col[1],$col[2],$col[3],$col[4],$col[5],$col[10],$col[15]\n";
		}
		else
		{
			$passedexons = $passedexons . "PASS,$col[0],$col[1],$col[2],$col[3],$col[4],$col[5],$col[10],$col[15]\n";
		}
	}
	
}

print "filter,gene,refSeq,exon,chr,startPos,endPos,min,pctGT$threshold\n";
print $failedexons;
print $passedexons;

close (COVERAGEFILE) or die ($!);
