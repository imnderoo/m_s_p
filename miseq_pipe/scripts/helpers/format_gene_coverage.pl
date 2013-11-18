#!/usr/bin/perl -w

use strict;
use warnings;
use File::Basename;

if ($#ARGV < 3)
{
	print "This tool parses through the coverage summary and extract the coverage stats for each gene in the gene list.\n";
	print "Usage format_coverage.pl <coverage_summary> <gene_list> <output_directory> <coverage depth threshold>\n";
	exit (1);
}

my $coverage = $ARGV[0];
my $geneListFile = $ARGV[1];
my $outDir = $ARGV[2];
my $threshold = $ARGV[3]; #for printing header mainly

my $covbase = basename($coverage);

$covbase =~  s/\.[^.]+$//;
$covbase =~ s/\./\_/g;

open (COVERAGEFILE, $coverage) or die($!);
open (GENELISTFILE, $geneListFile) or die($!);

mkdir "$outDir" unless (-d $outDir);

my %geneList = ();

while(my $line = <GENELISTFILE>) {
	chomp($line); #each line is a geneName
	my @genes = split('\t', $line); 
	$genes[0] =~ s/\r//g;
	$genes[0] =~ s/\s//g;
	$geneList{$genes[0]} = ();
}

while(my $line = <COVERAGEFILE>) {
	chomp($line);
	
	if ($line =~ m/(.+?)\t.+/) {
		my $geneName = $1;

		if (exists($geneList{$geneName})) {
			my @col = split("\t", $line);
			my $exon = $col[2];
			$exon =~ m/exon_([\d|\.]+)/;
			my $exonNum = $1;
			
			if(exists($geneList{$geneName}{$exonNum})) {
				$exonNum = $exonNum + 0.1;
			}

			$line =~ s/exon_([\d|\.]+)/exon_$exonNum/;

			$geneList{$geneName}{$exonNum} = "$line";
		}
	}	
}

foreach my $gene (keys %geneList) {

	open (OUT_FILE, ">$outDir/$covbase" . "_$gene.txt") or die $!;
	print OUT_FILE "gene\trefSeq\texon\tchr\tstartPos\tendPos\ttotal\tavg\tsample.total\tsample.avg\tmin\tq1\tmedian\tq3\tmax\tpctGT$threshold\tpctEQ0\n"; #Print header

	foreach my $exon (sort {$a <=> $b} keys %{$geneList{$gene}}) {
		my $line = $geneList{$gene}{$exon};
		print OUT_FILE "$line\n";
	}

	close (OUT_FILE);
}

close (COVERAGEFILE) or die ($!);
close (GENELISTFILE) or die ($!);
