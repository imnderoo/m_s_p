#!/usr/bin/python
# Scripts for converting .BED files to .interval_list files
# BED files are used by the coverage scripts while .interval_list are used by auto-classification

import os
import argparse
from bs4 import BeautifulSoup
from lxml import etree

def main():

	parser = argparse.ArgumentParser(description='Add Run info columns (SampleID, Flowcell, Run Date) to the front of the input files')
	parser.add_argument('runFolder', help='runFolder')
	parser.add_argument('inputFile', help='the temporary input file .tmp')

	args = parser.parse_args()
	
	runFolder = args.runFolder.rstrip('\/')
	inFileName = args.inputFile

	runName = os.path.basename(runFolder)

	# Parse run info from runFolder/Summary.xml
	summarySoup = createSummarySoup(runName,runFolder)
	runFolder = parseRunFolder(summarySoup)
	runSplit = runFolder.split("_")
	runDate = runSplit[0]
	runFlowCell = runSplit[3].split("-")[1]

	# Get sample id from CSV file
	fileID = os.path.basename(inFileName).split(".")[0]
	sampleID = fileID.split("_")[0]
	
	# Add run info to start of CSV file
	inCSV = open(inFileName, 'r')
	outCSV = open(inFileName.replace("tmp", "csv"), 'w')

	for i, line in enumerate(inCSV):
		if i == 0:
			line = "SampleID,Flowcell,RunDate," + line			
		else:
			line = sampleID + "," + runFlowCell + "," + runDate + "," + line

		outCSV.write(line)

	inCSV.close()
	outCSV.close()

def createSummarySoup(runName, runFolder):

	summaryFile = runFolder + "/Summary.xml"
	soup = BeautifulSoup(open(summaryFile), "xml")

	return soup

def parseRunFolder(soup):
	summary = soup.findAll("Summary")[0]
	chipSummary = summary.findAll("ChipSummary")[0]
	runFolder = chipSummary.findAll("RunFolder")[0].string

	return runFolder

main()

