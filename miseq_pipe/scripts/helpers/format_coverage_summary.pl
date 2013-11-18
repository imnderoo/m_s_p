#!/usr/bin/perl -w

use strict;
use warnings;

print "$#ARGV";

if ($#ARGV < 3)
{
	print "This tool takes the initial coverage summary file and adds annotations to it (Gene Name and Exon Name)\n";
	print "Usage format_coverage_summary.pl <Coverage_Interval_Summary> <RefGene> <RefSeq Exons> <Threshold>\n";
	exit;
}

#my $refGene = "/data1/genelists/refGene.broad.sorted.txt";
#my $refExon = "/data1/regions/refSeq.exon.sorted.bed"; # Or a selected interval list

my $coverage = $ARGV[0];
my $refGene = $ARGV[1];
my $refExon = $ARGV[2];
my $threshold = $ARGV[3];

open (REFGENEFILE, $refGene) or die($!);
open (REFEXONFILE, $refExon) or die($!);
open (COVERAGEFILE, $coverage) or die($!);
 
my %refGene = ();
my %refExon = ();

while(my $line = <REFGENEFILE>) {
	chomp($line);
	
	my @col = split(' ', $line);

	my $refSeqID = $col[1];
	my $geneName = $col[12];

	$refGene{"$refSeqID"}="$geneName";
}

while(my $line = <REFEXONFILE>) {
	chomp($line);
	#print "$line\n";

	my @col = split(' ', $line);

	my $chr = $col[0];
	my $startPos = $col[1]; # Need to add a +1 if using GATK
	my $endPos = $col[2];
	my $exonID = $col[3];
	
	my $chrPos = "$chr:$startPos-$endPos";

	my $refGeneID = "";
	my $exonName = "";

	if ($exonID =~ m/(NM|NR)_(.+?)_(.+?)_(.+?)_chr.+/)
	{
		$refGeneID = $1 . "_" . $2;
		$exonName = $3 . "_" . $4;
	}
	else
	{
		print "ERROR: $exonID\n";
	}

	my $geneName = "unknown";

	if (exists $refGene{"$refGeneID"}) {
		$geneName = $refGene{"$refGeneID"};
	}	

	# Instead of just using chrPos as the key, extract the exon number and add it to chrPos as such:
	# exonNumb.chrPos (treating it as a float during comparison)
	# so that the resulting graphs will be sorted in ascending order based on exon num first, then chr position
	
	$refExon{"$chrPos"} ="$geneName\t$refGeneID\t$exonName";
		
}

while(my $line = <COVERAGEFILE>) {
	chomp($line);
	
	if ($line =~ m/(chr.+?):(.+?)-(.+?)\t(.+)/) #$1 = chromosom $2 = startPos $3 = endPos $4 = coverageStats
	{
		if(exists $refExon{"$1:$2-$3"}) {
			my $exonInfo = $refExon{"$1:$2-$3"};
	        	$line = "$exonInfo\t$1\t$2\t$3\t$4"; #exonInfo contains: "GeneName  RefGeneID  ExonName" and is retrieved by chromosom position
		}
		else
		{
			$line = "unknown\tunknown\tunknown\t$1\t$2";
#			$unmatched = $unmatched + 1;
		}
	}
	else
	{
		$line = "gene\trefSeq\texon\tchr\tstartPos\tendPos\ttotal\tavg\tsample.total\tsample.avg\tmin\tq1\tmedian\tq3\tmax\tpctGT$threshold\tpctEQ0";
	}	

	print "$line\n";
}


close (REFEXONFILE) or die ($!);
close (COVERAGEFILE) or die ($!);
close (REFGENEFILE) or die ($!);
