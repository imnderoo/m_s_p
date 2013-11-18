#!/usr/bin/python
# Scripts for converting .BED files to .interval_list files
# BED files are used by the coverage scripts while .interval_list are used by auto-classification

import fileinput
import argparse
import vcreports.parser
import os
import geneTools.chr2cdsCoords as chr2cds

def main():
	parser = argparse.ArgumentParser(description='Flatten AnnoVar output (multianno.txt) such that LJB2_ALL, INFO, and GENO values are split into individual columns ')
	parser.add_argument('av', help='Path to annovar file')
	parser.add_argument('geneListFile', help='Path to genelist with accession #')
	parser.add_argument('variantsDBFile', help="Path to variantsDB")
	parser.add_argument('exonCoordsFile', help="Path to exonCoords")
	parser.add_argument('avFieldsList', help='AV Fields to extract')
	args = parser.parse_args()
	
	avIn = args.av
	avOut = os.path.splitext(args.av)[0] + "_filter.csv"
	geneList = getAccessionNumbers(args.geneListFile)
	variantsDict = getVariantsDB(args.variantsDBFile)
	avFieldsList = args.avFieldsList.split(',')
	exonDict = chr2cds.ucscGenesToDict(args.exonCoordsFile)

	filterAV(avIn, avOut, geneList, variantsDict, exonDict, avFieldsList)

def getAccessionNumbers(geneListFile):
	geneList = []
	infile = open(geneListFile)
	
	for line in infile:
		line = line.rstrip('\n\r')
		lineSplit = line.split('\t')
		geneList.append(lineSplit[1])

	return geneList

def getVariantsDB(variantsDBFile):
	variantsDict = {}
	infile = open (variantsDBFile)
	
	for line in infile:
		line = line.rstrip('\n\r')
		lineSplit = line.split('\t')
		# Use the original variant format as the variant key
		variantKey = lineSplit[0] + '\t' + lineSplit[2]
		variantsDict[variantKey] = '\t'.join(lineSplit[2:])

	return variantsDict

def filterAV(inFileName, outFileName, geneList, variantsDict, exonDict, avFieldsList):
	infile = open(inFileName)
	outfile = open(outFileName, 'w')

	# Keeps track of which columns to list (based on avFieldsList)
	colIdxDict = {}

	for i, line in enumerate(infile):
		line = line.rstrip('\n\r')
		lineSplit = line.split('\t')

		if i == 0:
			for j, heading in enumerate(lineSplit):
				if heading in avFieldsList:
					colIdxDict[j] = heading
		if i > 0:
			break
	
	infile.close()
	
	infile = open(inFileName)
	
	# Print body
	for i, line in enumerate(infile):
		line = line.rstrip('\n\r')
		lineSplit = line.split('\t')
		
		# Reset line to blank
		line = ""
		chrCoord = lineSplit[1]
		ref = lineSplit[3]
		alt = lineSplit[4]
		geneName = lineSplit[6].split('(')[0]
		geneName = geneName.split(';')[0]

		cDNA = ""
		cDNAflip = "" # This is necessary because some of the ref and alt in the variant DB is actually flipped

		for j in sorted(colIdxDict.keys()):

			heading = colIdxDict[j]
			lineSplit[j] = lineSplit[j].replace(",", ";")

			if heading == "Gene.refGene":
				if i == 0:
					line = line + "Gene.refGene" + ','
				else:
					line = line + geneName + ','
			elif heading == "AAChange.refGene":
				if i == 0:
					line = line + "AccessID.refGene,Exon.refGene,cDNAChange.refGene,AAChange.refGene" + ','
				else:

					if (exonDict[geneName]['strand'] == "+"):
						cdnaDict = chr2cds.find_cdna_fwd_strand_coords(exonDict, geneName, chrCoord)
					if (exonDict[geneName]['strand'] == "-"):
						cdnaDict = chr2cds.find_cdna_rev_strand_coords(exonDict, geneName, chrCoord)
						
					accession = cdnaDict["accession"]
					exon = cdnaDict["exon"]

					# Formatting allele for cDNAChange nomenclature
				
					# Reversre the nucleotide if the CDS is on the complementary strand
					if (exonDict[geneName]['strand'] == "-"):
					
						# Reference allele
						if ref == "G":
							ref = "C"
						elif ref == "C":
							ref = "G"
						elif ref == "A":
							ref = "T"
						elif ref == "T":
							ref = "A"

						# Alternate allele
						if alt == "G":
							alt = "C"
						elif alt == "C":
							alt = "G"
						elif alt == "A":
							alt = "T"
						elif alt == "T":
							alt = "A" 

					if ref == "-":
						cDNA = cdnaDict["cdna"] + "ins" + alt
					elif alt == "-":
						cDNA = cdnaDict["cdna"] + "del" + ref
					else:
						cDNA = cdnaDict["cdna"] + ref + ">" + alt
						cDNAflip = cdnaDict["cdna"] + alt + ">" + ref
					
					line = line + ','.join([accession, exon, cDNA, "NA"]) + ","

					"""
					if lineplit[j] == "NA" or lineSplit[j] == "UNKNOWN":
						line = line + ','.join(["NA","NA","NA","NA"]) + ','
						cDNA = "NA"

					else:
						transcripts = lineSplit[j].split(',')
						transcriptFound = 0;

						for transcript in transcripts:
							values = transcript.split(':')
							cDNA = values[3]							

							if values[1] in geneList:
								values.pop(0)
								line = line + ','.join(values) + ','
								transcriptFound = 1;

						if transcriptFound == 0:
								line = line + ','.join(["NA","NA","NA","NA"]) + ','
						"""

			elif heading == "AD":
				if i == 0:
					line = line + "AD.REF,AD.ALT" + ','
				else:
					adSplit = lineSplit[j].split(';')
				
					# Split the AD intp REF and ALT and then add it.
					line = line + adSplit[0] + ',' + adSplit[1] + ','
			else:
				line = line + lineSplit[j] + ','

		# Update Variant Fields Based on Variant DB
		if i == 0:
			line = line + "cDNAChange,AAChange,cDNAType,AAType,LatestCategory,ACMG,LatestDate,ReportCat" + ","
		else:
			variantKey = geneName + "\t" + cDNA
			variantKeyflip = geneName + "\t" + cDNAflip # Used to chek in case the ref and alt is flipped in vDB
			acmg = 0

			# Update Variant Fields from Variant DB
			if variantKey in variantsDict:
				acmg = variantsDict[variantKey].split("\t")[5]
				cdnaType = variantsDict[variantKey].split("\t")[2]
				aaType = variantsDict[variantKey].split("\t")[3]
				
				if acmg == '':
					acmg = 0
					line = line + "-,-,-,-,Not Assessed,-,-" + ","
				else:
					line = line + variantsDict[variantKey].replace("\t",",") + ","
				
			elif variantKeyflip in variantsDict:
				acmg = variantsDict[variantKeyflip].split("\t")[5]
				cdnaType = variantsDict[variantKeyflip].split("\t")[2]
				aaType = variantsDict[variantKeyflip].split("\t")[3]
			
				if acmg == '':
					acmg = 0
					line = line + "-,-,-,-,Not Assessed,-,-" + ","
				else:
					line = line + variantsDict[variantKeyflip].replace("\t",",") + ","

			else:
				line = line + "-,-,-,-,Not Assessed,-,-" + ","
			
			# Update Report Cat for Jasper Reports
			if int(acmg) == 0:
				line = line + "1" + "," # 1 - Not Assessed
			elif int(acmg) <= 3:
				line = line + "0" + "," # 0 - Confirm 
			elif int(acmg) == 4:
				line = line + "2" + "," # 2 - Prediced Benign
			else:
				line = line + "3" + "," # 3 - Benign
		
		line = line + '\n'
		line = line.replace(',\n','\n')
			
		outfile.write(line)

	infile.close()
	outfile.close()

main()
