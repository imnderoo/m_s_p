#!/usr/bin/perl -w

use strict;
use warnings;
use POSIX;

if ($#ARGV < 1)
{
	print "This tool takes the histogrm and outputs summary statistics such as min, q1, median, q3, max, pct>threshold. The output is sent to STDPUT.\n";
	print "Usage create_coverage_summary.pl <coverage.hist> <coverage depth threshold>\n";
	exit; 
}

my $coverage = $ARGV[0];
my $threshold = $ARGV[1];
#my $threshold = 20;

open (COVERAGEFILE, "$coverage") or die($!);
 
my %interval = ();

my %results = ();
my @depthArray = ();
my $prevExon = "empty.0.0.exon"; #The funny name is to follow the format of exonKey "chr.startPos.endPos.exonName"
my $total = 0;
my $numAboveThreshold = 0;
my $numAtZero = 0;

while(my $line = <COVERAGEFILE>) {
	chomp($line);
	
	my ($chr,$startPos,$endPos,$exonName,$depth,$numAtDepth,$exonKey);

	if( $line =~ /^chr.+/) {
		
		my @col = split(' ', $line);

		$chr= $col[0];
		$startPos = $col[1];
		$endPos = $col[2];
		$exonName = $col[3];
		$exonName =~ m/.+?_exon_(.+?)_chr.+/;

		$depth = $col[6];
		$numAtDepth = $col[7];

		$exonKey =  "$chr.$startPos.$endPos.$exonName"; 	
	}

	if ( $line =~ /^all.+/) {
		my @col = split(' ', $line);

		$chr= "all";
		$startPos = 0;
		$endPos = 0;
		$exonName = "exon";

		$depth = $col[1];
		$numAtDepth = $col[2];

		$exonKey =  "$chr.$startPos.$endPos.$exonName"; 	
	}
	
	if(!exists($interval{$exonKey})) {

   	  	my @prevPosCol = split ('\.', $prevExon);
   	    my $prevChrPos = $prevPosCol[0];
		my $prevStart = $prevPosCol[1];
		my $prevEnd = $prevPosCol[2];
		my $prevExonNum = $prevPosCol[3];

		# Used to remove deplicate exon entries (This happens when a gene have multiple transcripts. Each transcript shows u as a separate exon)
		
		if(!exists($results{$prevChrPos}{$prevStart})) {

         	my $numBase = $#depthArray + 1;
				
			if($prevExon eq "empty.0.0.exon") {
				print "chr\ttotal\tavg\tsample_total\tsample_avg\tmin\tq1\tmedian\tq3\tmax\tpctGT$threshold\tpctEQ0\n";
				$results{"empty"}{0}="";
			}
			elsif($prevExon eq "all.0.0.exon") {
				# Do nothing
			}
			else {
            	# Calculate Q1, Median, Q3. 
                # Uses floor function because fraction of a base doesn't make sense. Uses floor b/c the depth is for QC purposes - better safe than sorry.
                # The -1 is for array index
		        my $q1=floor(1 / 4 * ($numBase + 1)) - 1;
                my $q2=floor(1 / 2 * ($numBase + 1)) - 1;
                my $q3=floor(3 / 4 * ($numBase + 1)) - 1;
                my $avg=floor($total / ($numBase));
				my $pctAboveThreshold= sprintf ("%.2f", $numAboveThreshold / $numBase * 100);
				my $pctZero = sprintf ("%.2f", $numAtZero / $numBase * 100);
				# These two variables are placeholders in order to be consistent with GATK DepthOfCoverage output
				my $sampleTotal = $total;
				my $sampleAvg = $avg;

				my @sortedArray = sort {$a <=> $b} @depthArray;			

				$results{$prevChrPos}{$prevStart} = "$prevChrPos:$prevStart-$prevEnd\t$total\t$avg\t$sampleTotal\t$sampleAvg\t$sortedArray[0]\t$sortedArray[$q1]\t$sortedArray[$q2]\t$sortedArray[$q3]\t$sortedArray[$#sortedArray]\t$pctAboveThreshold\t$pctZero\n";
			}
		}
	
		$interval{$exonKey} = "";
	
		# Reset variables
		@depthArray = ();
		$numAtZero = 0;
		$numAboveThreshold = 0;
		$total = 0;
	}	
	

	$prevExon = $exonKey;

	# For each number of base, add depth to the beginning of array.

	for (my $count = 1; $count <= $numAtDepth ; $count++) {
                push(@depthArray, $depth);
	}

	$total = $total + ($depth * $numAtDepth);
	
	if($depth == 0) {
		$numAtZero = $numAtZero + $numAtDepth;
	}	

	if($depth >= $threshold) {
		$numAboveThreshold = $numAboveThreshold + $numAtDepth;
	} # end line check	
} # end While

foreach my $chr_key(sort keys %results) {
	foreach my $pos_key(sort { $a <=> $b} keys %{$results{$chr_key}}) {
 		print "$results{$chr_key}{$pos_key}";
	}
}

close (COVERAGEFILE) or die ($!);
