#!/usr/bin/python
# Scripts for converting .BED files to .interval_list files
# BED files are used by the coverage scripts while .interval_list are used by auto-classification

import os, time, sys, gc
import re
import argparse
import csv
import bisect

def ucscGenesToDict(exonFileName):

	# Parse exonic info
	csvReader = csv.DictReader(open(exonFileName, 'r'), delimiter="\t")
	
	genesDict = {}

	for row in csvReader:
		geneName = row["hg19.refGene.name2"]
		accession = row["hg19.refGene.name"]
		exonStarts = row["hg19.refGene.exonStarts"]
		exonEnds = row["hg19.refGene.exonEnds"]
		strand = row["hg19.refGene.strand"]
		cds = row["hg19.cds.name"]

		genesDict[geneName] = {}
		genesDict[geneName]['accession'] = accession

		exonStartList = exonStarts.split(',')
		exonEndList = exonEnds.split(',')
		
		cdsStart = cds.split("..")[0]
		cdsEnd = cds.split("..")[1]
		
		exonStartList.pop()
		exonEndList.pop()

		#genesDict[geneName]['exonStartList'] = sorted(map(int, exonStartList.pop()))
		#genesDict[geneName]['exonEndList'] = sorted(map(int, exonEndList.pop()))
		
		genesDict[geneName]['exonStartList'] = sorted(map(int, exonStartList))
		genesDict[geneName]['exonEndList'] = sorted(map(int, exonEndList))
		genesDict[geneName]['strand'] = strand
		genesDict[geneName]['cdsStart'] = int(cdsStart)
		genesDict[geneName]['cdsEnd'] = int(cdsEnd)

	return genesDict

def find_fwd_strand_exon_start_idx(exonStartList, coord):
	'Find rightmost exonStartCoord index that is less than x'
	idx = bisect.bisect_right(exonStartList, int(coord))
	
	if idx:
		return idx-1
	else:
		return -1
		
def find_rev_strand_exon_start_idx(exonStartList, coord):
	'Find leftmost exonStartCoord index that is greater than or equal to x'
	idx = bisect.bisect_left(exonStartList, int(coord))
	
	if idx != len(exonStartList):
		return idx
	else:
		return -1

def find_cdna_rev_strand_coords(genesDict, geneName, coord):
	cdnaDict = {}
	coord = int(coord)
	cdnaDict["accession"] = genesDict[geneName]["accession"]
	cdnaDict["exon"] = "intron"
	cdnaDict["region"] = "unknown"
	cdnaDict["cdna"] = "unknown"

	# Flip the definition of Start and End when it is on the reverse strand
	exonStartList = genesDict[geneName]['exonEndList']
	# +1 because the Start coordinates are 0-based but the chromosome coordinates will be 1-based
	exonEndList = genesDict[geneName]['exonStartList']
	exonEndList = sorted(map(lambda x:x+1, exonEndList))
	#coord = int(args.coord)

	# Have to subtract by sumCDSCoord cdsStartCoords when displaying cDNA coordinates.
	# cDNA coordinate doesn't start with the first base of exon, but somewhere down the road.
	# cDNA start and end can be grabbed off of ucsc hg19.cds table, name field.
	cdsStartCoord = genesDict[geneName]['cdsStart']
	# cdsStartCoord = 0
	cdsEndCoord = genesDict[geneName]['cdsEnd']
	
	exonIdx = find_rev_strand_exon_start_idx(exonStartList, coord)

	if(exonIdx == -1):
		exonCoordsDiff = coord - exonStartList[0]
		cdnaDict["region"] = "Before transcription begin"
		cdnaDict["cdna"] = "c.-" + str(exonCoordsDiff)

	# Find the CDNA coord of the SNP Coord
	else:
		sumCDSCoord = 1 # First base of exon 1 starts at 1
		finalCDSCoord = 0

		# Find the CDS of the closest exonStart that's less than SNP Coord
		# Convert all the exon coordinates to CDS coordinates 
		for idx in range(len(exonStartList)-1, exonIdx, -1):
			# prevExonEndCDSCoord = (exonEndList[idx-1] - exonStartList[idx-1] + sumCDSCoord)
			# curExonStartCDSCoord = prevExonEndCDSCoord + 1
			sumCDSCoord = (exonStartList[idx] - exonEndList[idx]) + sumCDSCoord + 1
			
			#print (str(exonEndList[idx]) + "-" + str(exonStartList[idx]) + " " + str(sumCDSCoord))


		# Check if the coords is less (bc it's rev strand) than the exon end coords.
		# If it is, the SNP is intronic / splicing
		if (coord < exonEndList[exonIdx]):
			# If the SNP position is AFTER the last exon coordinate
			if (exonIdx == 0):
				exonCoordsDiff = exonEndList[exonIdx] - coord
				cdnaDict["region"] = "After transcription ends"
				cdnaDict["cdna"] = "c.*" + str(exonCoordsDiff)

			else:
				distFromExonEnd = exonEndList[exonIdx] - coord
				distFromNextExonStart = coord - exonStartList[exonIdx-1]
			
				# If it's closer to exon End, get CDS of exon end and + the difference (ex. c.134+3)
				if (distFromExonEnd < distFromNextExonStart):
					# Update sumCDS b/c CDNA is relative to cur. exon end CDS coords instead of cur.  exon start
					sumCDSCoord = (exonStartList[exonIdx] - exonEndList[exonIdx]) + sumCDSCoord
					finalCDSCoord = sumCDSCoord - cdsStartCoord + 1
					cdnaDict["region"] = "Intronic: Closer to exon end"
					cdnaDict["cdna"] = "c." + str(finalCDSCoord)  + "+" + str(distFromExonEnd)

				# If it's closer to the next exon Start, get the CDS of exon and - the difference (ex. c.134-4)
				else:
					# Update sumCDS b/c CDNA is relative to next exon start CDS coords instead of prev. exon start
					sumCDSCoord = (exonStartList[exonIdx] - exonEndList[exonIdx]) + sumCDSCoord + 1
					finalCDSCoord = sumCDSCoord - cdsStartCoord + 1
					cdnaDict["region"] = "Intronic: Closer to exon start."
					cdnaDict["cdna"] = "c." + str(finalCDSCoord) + "-" + str(distFromNextExonStart)
			
		# If it is less than the exon end coords, the SNP is exonic.
		else:
			finalCDSCoord = sumCDSCoord - cdsStartCoord + 1
			exonCoordsDiff = exonStartList[exonIdx] - coord

			cdnaDict["exon"] = str(len(exonEndList) - (exonIdx))
			cdnaDict["region"] = "Exonic"
			cdnaDict["cdna"] = "c." + str(finalCDSCoord + exonCoordsDiff)

	return cdnaDict

def find_cdna_fwd_strand_coords(genesDict, geneName, coord):

	cdnaDict = {}
	coord=int(coord)
	cdnaDict["accession"] = genesDict[geneName]["accession"]
	cdnaDict["exon"] = "-"
	cdnaDict["cdna"] = "unknown"
	cdnaDict["region"] = "unknown"
	
	exonStartList = genesDict[geneName]['exonStartList']
	# Add 1 to the exonStart because it's a 0-based coordinate system, but chromosome coordinate is 1-based
	exonStartList = sorted(map(lambda x:x+1, exonStartList))
	exonEndList = genesDict[geneName]['exonEndList']

	# Have to subtract by sumCDSCoord cdsStartCoords when displaying cDNA coordinates.
	# cDNA coordinate doesn't start with the first base of exon, but somewhere down the road.
	# cDNA start and end can be grabbed off of ucsc hg19.cds table, name field.
	cdsStartCoord = genesDict[geneName]['cdsStart']
	cdsEndCoord = genesDict[geneName]['cdsEnd']
	
	exonIdx = find_fwd_strand_exon_start_idx(exonStartList, coord)

	if(exonIdx == -1):
		exonCoordsDiff = exonStartList[0] - coord
		cdnaDict["region"] = "Before transcription starts"
		cdnaDict["cdna"] = "c.-" + str(exonCoordsDiff)

	# Find the CDNA coord of the SNP Coord
	else:
		sumCDSCoord = 1 # First base of exon 1 starts at 1

		finalCDSCoord = 0

		# Find the CDS of the closest exonStart that's less than SNP Coord
		# Convert all the exon coordinates to CDS coordinates 
		for idx in range(0, exonIdx):
			# prevExonEndCDSCoord = (exonEndList[idx-1] - exonStartList[idx-1] + sumCDSCoord)
			# curExonStartCDSCoord = prevExonEndCDSCoord + 1
			sumCDSCoord = (exonEndList[idx]-exonStartList[idx]) + sumCDSCoord + 1
			
			# print (str(exonEndList[idx]) + "-" + str(exonStartList[idx]) + " " + str(sumCDSCoord))

		# Check if the coords is greater than the exon end coords.
		# If it is, the SNP is intronic / splicing
		if (coord > exonEndList[exonIdx]):
			# If the SNP position is AFTER the last exon coordinate
			if (exonIdx == len(exonEndList)-1):
				exonCoordsDiff = coord - exonEndList[exonIdx]
				cdnaDict["region"] = "After transcription ends"
				cdnaDict["cdna"] = "c.*" + str(exonCoordsDiff)
			else:
				distFromExonEnd = coord - exonEndList[exonIdx]
				distFromNextExonStart = exonStartList[exonIdx+1] - coord
			
				# If it's closer to exon End, get CDS of exon end and + the difference (ex. c.134+3)
				if (distFromExonEnd < distFromNextExonStart):
					# Update sumCDS b/c CDNA is relative to cur. exon end CDS coords instead of cur.  exon start
					sumCDSCoord = (exonEndList[exonIdx] - exonStartList[exonIdx]) + sumCDSCoord
					finalCDSCoord = sumCDSCoord - cdsStartCoord + 1
					cdnaDict["region"] = "Intronic: Closer to exon end"
					cdnaDict["cdna"] = "c." + str(finalCDSCoord)  + "+" + str(distFromExonEnd)

				# If it's closer to the next exon Start, get the CDS of exon and - the difference (ex. c.134-4)
				else:
					# Update sumCDS b/c CDNA is relative to next exon start CDS coords instead of prev. exon start
					sumCDSCoord = (exonEndList[exonIdx] - exonStartList[exonIdx]) + sumCDSCoord + 1
					finalCDSCoord = sumCDSCoord - cdsStartCoord + 1
					cdnaDict["region"] = "Intronic: Closer to exon start"
					cdnaDict["cdna"] = "c." + str(finalCDSCoord) + "-" + str(distFromNextExonStart)

		# If it is less than the exon end coords, the SNP is exonic.
		else:
			finalCDSCoord = sumCDSCoord - cdsStartCoord + 1
			exonCoordsDiff = coord - exonStartList[exonIdx]

			cdnaDict["region"] = "Exonic"
			cdnaDict["cdna"] = "c." + str(finalCDSCoord + exonCoordsDiff)
			cdnaDict["exon"] = str(exonIdx + 1)

			# TODO: If the finalCDSCoord is negative, then treat it as BEFORE transcription.
			# TODO: Before transcription means taking the transcription start

	return cdnaDict

