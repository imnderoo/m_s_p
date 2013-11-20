#!/usr/bin/python
# Scripts for converting .BED files to .interval_list files
# BED files are used by the coverage scripts while .interval_list are used by auto-classification

import argparse
import csv

def formatVarSummary(csvReader):
	# Before sorting, remove variants that have not passed QC
	# Also, add a blank dummy line to group categories (defined for jasperReport) so that it will be reported as "None Detected"
	csvResult = []
	reportCatCountDict = {'0':0, '1':0, '2':0, '3':0}
	# 0 - Confirm
	# 1 - Not Assessed
	# 2 - Predicted Benign
	# 3 - Benign

	dummyLine = {}

	for row in csvReader:
		# Dummy line for later
		for keys in row:
			dummyLine[keys] = row[keys]

		# Exclude variants that did not pass filter
		qc = row["FILTER"]
		if qc == "PASS":
			csvResult.append(row)
			# Track number of passed variants in each category
			reportCat = row["ReportCat"] 
			reportCatCountDict[reportCat] = reportCatCountDict[reportCat] + 1

	for reportCatKey in reportCatCountDict:
		if reportCatCountDict[reportCatKey] == 0:
			dummyLineInstance = {}
			for dummyKey in dummyLine:
				dummyLineInstance[dummyKey] = "-"
			
			dummyLineInstance["ReportCat"] = reportCatKey
			csvResult.append(dummyLineInstance)

	return csvResult

def formatExonSummary(csvReader):
	# Add a blank dummy line to group categories (defined for jasperReport) so that it will be reported as "None Detected"
	csvResult = []
	reseqCatCountDict = {'PASS':0, 'FAIL':0}
	dummyLine = {}

	for row in csvReader:
		# Dummy line for later
		for keys in row:
			dummyLine[keys] = row[keys]

		# Exclude variants that did not pass filter
		csvResult.append(row)
		# Track number of passed variants in each category
		reseqCat = row["filter"] 
		reseqCatCountDict[reseqCat] = reseqCatCountDict[reseqCat] + 1

	for reseqCatKey in reseqCatCountDict:
		if reseqCatCountDict[reseqCatKey] == 0:
			dummyLineInstance = {}
			for dummyKey in dummyLine:
				dummyLineInstance[dummyKey] = "-"
			
			dummyLineInstance["filter"] = reseqCatKey
			csvResult.append(dummyLineInstance)

	return csvResult

def main():

	parser = argparse.ArgumentParser(description='Sort the input file according to the input list of column headers ')
	parser.add_argument('inputFile', help='the temporary input file')
	parser.add_argument('sortFields', help='fields to sort. Separate by commas. List in order of sort')

	args = parser.parse_args()
	
	inFileName = args.inputFile
	
	inCSV = open(inFileName, 'r')
	csvReader = csv.DictReader(open(inFileName, 'r'))
	
	sortFields = args.sortFields
	sortFields = sortFields.split(",")[::-1]

	csvResult = csvReader

	# Add a special case for Variant Summary and CSV formatting.
	# 1. Filter out variants that did not pass filter
	# 2. Add a dummy row if there are 0 variants in a specific category. \
	# This is a workaround needed for JasperReport to print empty categories

	if "varsummary" in inFileName:
		csvResult = formatVarSummary(csvReader)
	if "exonsummary" in inFileName:
		csvResult = formatExonSummary(csvReader)

	for fields in sortFields:
		csvResult = sorted(csvResult, key=lambda d: d[fields])		

	inFileExt="." + inFileName.split(".")[-1]
	newInFileExt="_sorted" + inFileExt

	csvWriter = csv.DictWriter(open(inFileName.replace(inFileExt, newInFileExt), 'w'), csvReader.fieldnames)
	csvWriter.writeheader()
	csvWriter.writerows(csvResult)

	inCSV.close()

main()

