#!/usr/bin/python

import os
import argparse
from decimal import *
from bs4 import BeautifulSoup
from lxml import etree
from illuminate import InteropDataset

# Mainly need to update the path to DATA/SUMMARY as well as the path to BAM and VCF files.
# InterOp, SampleSheet, and RunInfo is still access from RunFolder

def parseSampleSheet(sampleSheetName, sampleID):
	sampleSheetDict = {}

	sampleSheet = open(sampleSheetName, 'r')
	lines = sampleSheet.readlines()

	manifestDict = {}

	for line in enumerate(lines):
		if "Investigator" in line[1]:
			sampleSheetDict['investigator'] = line[1].split(',')[1].rstrip('\r\n')
		elif "Assay" in line[1]:
			sampleSheetDict['assay'] = line[1].split(',')[1].rstrip('\r\n')
		elif "Experiment" in line[1]:
			sampleSheetDict['experiment'] = line[1].split(',')[1].rstrip('\r\n')
		elif "Description" in line[1]:
			sampleSheetDict['description'] = line[1].split(',')[1].rstrip('\r\n')
		elif "[Manifests]" in line[1]:
			i = 1
			manifestLine = lines[line[0]+i]
			while manifestLine not in ['\n', '\r\n']:
				manifestID = manifestLine.split(',')[0]
				manifestFile = manifestLine.split(',')[1].rstrip('\r\n')
				manifestDict[manifestID] = manifestFile
				i=i+1
				manifestLine = lines[line[0]+i]

		elif sampleID == line[1].split(',')[0]:
			sampleSheetDict['indexi5'] = line[1].split(',')[6] + ":" + line[1].split(',')[7]
			sampleSheetDict['indexi7'] = line[1].split(',')[4] + ":" + line[1].split(',')[5]
			manifestID = line[1].split(',')[8]
			sampleSheetDict['manifest'] = manifestDict[manifestID]

	return sampleSheetDict

def main():

	parser = argparse.ArgumentParser(description='Extract QC metric from MiSeq Summary.xml')
	parser.add_argument('runFolder', help='runFolder')
	parser.add_argument('inFileName', help='input file')
	parser.add_argument('variantdb', help='variantdb')
	parser.add_argument('script', help='script')

	args = parser.parse_args()
	
	runFolder = args.runFolder.rstrip('\/')
	runName = os.path.basename(runFolder)

	inFileName = args.inFileName

	# Get sample id from CSV file
	fileID = os.path.basename(inFileName).split(".")[0]
	sampleID = "_".join(fileID.split('_')[0:-1])
	sheetSampleID = fileID.split('_')[0]

	# Parse QC info from .xml
	summarySoup = createSummarySoup(runName,runFolder)

	baseYield = parseYield(summarySoup)
	run = parseRunFolder(summarySoup)

	# print (baseYield)

	runSplit = run.split("_")
	runDate = runSplit[0]
	runFlowCell = runSplit[3].split("-")[1]

	# print (runDate + " " + runFlowCell)

	outDir = os.path.dirname(inFileName)
	outFileName =  outDir + "/" + sampleID + "_runqc_sorted.csv"
	outCSV = open(outFileName, 'w')

	# Use Illuminate:InterOp to get stats
	
	myQC = InteropDataset(runFolder)

	tileQC = myQC.TileMetrics()
	density = str(tileQC.mean_cluster_density_pf)
	pctClusterPF = str(tileQC.percent_pf_clusters)

	qualQC = myQC.QualityMetrics()
	q30List = qualQC.read_qscore_results['q30']
	pctGT30 = str(sum(q30List) / len(q30List))

	indexQC = myQC.IndexMetrics()
	readsPF = str(indexQC.total_ix_reads_pf)
	
	# Parse sample sheet
	sampleSheet = runFolder + "/SampleSheet.csv"
	sampleSheetDict = {}
	sampleSheetDict = parseSampleSheet(sampleSheet, sheetSampleID)

	# Print Header

	header="SampleID,Flowcell,RunDate,yieldTotal(G),pctGT30,density(k/mm2),pctClusterPF,readsPF(M),script,manifest,variantdb,indexi5,indexi7,investigator,description,assay,experiment"
	
	manifest = sampleSheetDict['manifest']
	indexi5 = sampleSheetDict['indexi5']
	indexi7 = sampleSheetDict['indexi7']
	investigator = sampleSheetDict['investigator']
	description = sampleSheetDict['description']
	assay = sampleSheetDict['assay']
	experiment = sampleSheetDict['experiment']

	metrics = ",".join([sheetSampleID,runFlowCell,runDate,baseYield,pctGT30,density,pctClusterPF,readsPF,args.script,manifest,args.variantdb,indexi5,indexi7,investigator,description,assay,experiment])
	
	outCSV.write(header + "\n")
	outCSV.write(metrics + "\n")

	outCSV.close()

def createSummarySoup(runName, runFolder):

	summaryFile = runFolder + "/Data/Intensities/BaseCalls/Alignment/Summary.xml"
	soup = BeautifulSoup(open(summaryFile), "xml")

	return soup

def parseYield(soup):

	summary = soup.findAll("Summary")[0]
	
	chipResultsSummary = summary.findAll("ChipResultsSummary")[0]

	baseYield = chipResultsSummary.findAll("yield")[0].string
	
	return baseYield

def parseRunFolder(soup):
	summary = soup.findAll("Summary")[0]
	chipSummary = summary.findAll("ChipSummary")[0]
	runFolder = chipSummary.findAll("RunFolder")[0].string

	return runFolder

main()

