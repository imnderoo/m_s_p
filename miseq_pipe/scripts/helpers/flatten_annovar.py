#!/usr/bin/python
# Scripts for converting .BED files to .interval_list files
# BED files are used by the coverage scripts while .interval_list are used by auto-classification

import fileinput
import argparse
import vcreports.parser
import os
def main():
	parser = argparse.ArgumentParser(description='Flatten AnnoVar output (multianno.txt) such that LJB2_ALL, INFO, and GENO values are split into individual columns ')
	parser.add_argument('av', help='Path to Annovar file')
	parser.add_argument('infoFields', help='INFO Fields to extract')
	parser.add_argument('genoFields', help='GENO Fields to extract')

	args = parser.parse_args()
	
	avIn = args.av
	avOut = os.path.splitext(args.av)[0] + "_flatten.txt"

	infoFieldsList = args.infoFields.split(',')
	genoFieldsList = args.genoFields.split(',')
	
	flattenAV(avIn, avOut, infoFieldsList, genoFieldsList)

def flattenAV(inFileName, outFileName, infoFieldsList, genoFieldsList):
	infile = open(inFileName)
	outfile = open(outFileName, 'w')

	# Quickly jump to line 2 to read the file headers
	header = ""

	for i, line in enumerate(infile):
		
		line = line.rstrip('\n\r')
		lineSplit = line.split('\t')		

		if i == 0:
			header = line

			# Replace if these headers are kept
			header = header.replace("ljb2_all", "SIFT\tPP2-HVar\tPP2-HVar-Category\tPP2-HDiv\tPP2-HDiv-Category\tLRT\tLRT-Category\tMT\tMT-Category\tMA\tMA-Category\tFATHMM\tGERP++\tPhyloP\tSiPhy")
			header = header.replace("Otherinfo", "Zygosity\tCHR\tCHRPOS\tdbSNPID\tREF\tOBS\tQUAL\tFILTER")
		if i == 1:
			infoCol = lineSplit[-3]
			genoHeader = lineSplit [-2]
			genoCol = lineSplit[-1]

			infoDict = vcreports.parser.parseAndFilterInfoValues(infoCol, infoFieldsList)
			genoDict = vcreports.parser.parseAndFilterGenoValues(genoHeader, genoCol, genoFieldsList)	
			
			for infoKey in sorted(infoDict.keys()):
				header = header + "\t" + infoKey

			for genoKey in sorted(genoDict.keys()):
				header = header + "\t" + genoKey

		if i > 1:
			break
	
	infile.close()
	
	outfile.write(header + "\n")	
	
	infile = open(inFileName)
	
	# Print body
	for i, line in enumerate(infile):
		if i != 0:
			line = line.rstrip('\n\r')
			lineSplit = line.split('\t')
		
			# Join everything that doesn't need post-processing

			
			line = '\t'.join(lineSplit[0:-12])
			
			# Split LJB Scores into individual columns
			#ljbChunk = lineSplit[-12]
			#ljbCols = ""
			#if (ljbChunk == "NA"):
			#	ljbCols = ".\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t"
			#else:
			#	ljbCols = ljbChunk.replace(",","\t")
			#	ljbCols = ljbCols + '\t'

			# Join the Misc info columns
			miscCols = '\t'.join(lineSplit[-12:-3])
			miscCols = miscCols + '\t'

			# Split INFO column into individual fields
			infoChunk = lineSplit[-3]
			infoCols = ""
			infoDict = vcreports.parser.parseAndFilterInfoValues(infoChunk, infoFieldsList)
			
			for infoKey in sorted(infoDict.keys()):
				infoCols = infoCols + infoDict[infoKey] + "\t"
			
			# Split GENO Columns into individual fields
			genoHeadChunk = lineSplit[-2]
			genoValChunk = lineSplit[-1]
			genoCols = ""
				
			genoDict = vcreports.parser.parseAndFilterGenoValues(genoHeadChunk, genoValChunk, genoFieldsList)

			for genoKey in sorted(genoDict.keys()):
				genoCols = genoCols + genoDict[genoKey] + "\t"
				
			genoCols = genoCols.replace("/",",") # This fixes the GT field (1/1) which becomes a date when uploaded to excel
			
			line = line + "\t" + miscCols + infoCols + genoCols + "\n"

			line = line.replace("\t\n","\n")
			
			outfile.write(line)

	infile.close()
	outfile.close()

main()
