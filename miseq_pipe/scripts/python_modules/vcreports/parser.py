#!/usr/bin/python

import re

def parseGeneList(genelistFile):
	infile = open(genelistFile)

	genes = []

	for line in infile:
		geneName = line.rstrip('\n\r')
		geneName = geneName.strip()
		genes.append(geneName)

	infile.close()
	return genes

# This method parses the MiSeqAV (AnnoVar) output
# The input geneList is used to filter out AnnoVar entries that are not part of the interested genes
# Returns a dictionary of dictionary: avDict['snp_pos_and_allele'] = avFieldDict
# where avFieldDict is a dictionary where: avFieldDict['field_name'] = field_value

def parseMasterVariants(masterVariantsFile):
	infile = open(masterVariantsFile)
	masterDict = {}
	#interestedFields = ["autoclass", "mutationType"]

	for line in infile:
		masterFieldDict = {}
		lineSplit = line.split('\t') #Variant Position
		chr = lineSplit[0].replace('chr', '').replace('"', '')
		chrPos = lineSplit[1].replace('"', '')
		ref = lineSplit[2].replace('"', '')
		alt = lineSplit[3].replace('"', '')
		
		variantPos = chr + "\t" + chrPos + "\t" + ref + "\t" + alt 

		#masterFieldDict["rsID"] = lineSplit[4].replace('"', '')
		masterFieldDict["region"] = lineSplit[7].replace('"', '')
		masterFieldDict["gene"] = lineSplit[8].replace('"', '')
		#masterFieldDict["accession"] = lineSplit[9].replace('"', '')
		#masterFieldDict["cdnaChange"] = lineSplit[11].replace('"', '')
		#masterFieldDict["aaChange"] = lineSplit[12].replace('"', '')
		masterFieldDict["mutationType"] = lineSplit[14].replace('"', '')
		masterFieldDict["autoclass"] = lineSplit[20].replace('"', '')

		masterDict[variantPos] = masterFieldDict

	infile.close()
	return masterDict

# This method parses the MiSeqAV (AnnoVar) output
# The input geneList is used to filter out AnnoVar entries that are not part of the interested genes
# Returns a dictionary of dictionary: avDict['snp_pos_and_allele'] = avFieldDict
# where avFieldDict is a dictionary where: avFieldDict['field_name'] = field_value

def parseAndFilterMiSeqAV(avFile, interestedGeneList):
	infile = open(avFile)
	avDict = {}

	for line in infile:
	
		## Note that SB can appear in either Info or Geno fields. Therefore, it's dynamically added based on where it's found
		interestedInfoFields = ["DP", "MQ"]
		interestedGenoFields = ["GT", "GQ", "AD"]

		if not line.startswith("Func"):
			line = line.rstrip('\n\r')

			quoteCommas = re.findall('"[\:\/|\-|\.|\w|\=|\;]+,.+?"', line)

			if len(quoteCommas) > 0:
				for quoteComma in quoteCommas:
					fixedQuoteComma = quoteComma.replace(',', ' & ')
					line = line.replace(quoteComma, fixedQuoteComma)

			lineSplit = line.split(',')
			
			avFieldDict = {}
			avFieldDict["func"] = lineSplit[0].replace('"', '')
			avFieldDict["gene"] = lineSplit[1].replace('"', '')
			avFieldDict["exonicfunc"] = lineSplit[2].replace('"', '')
			avFieldDict["aaChange"] = lineSplit[3].replace('"', '')
			avFieldDict["popESP"] = lineSplit[6].replace('"', '')
			avFieldDict["pop1000g"] = lineSplit[7].replace('"', '')
			avFieldDict["dbsnpID"] = lineSplit[8].replace('"', '')
			avFieldDict["zygosity"] = lineSplit[26].replace('"', '')
			
			chr = lineSplit[21].replace("chr", "").replace('"', '')
			chrPos = lineSplit[22].replace('"', '')
			ref = lineSplit[24].replace('"', '')
			alt = lineSplit[25].replace('"', '')
				
			variantPos = chr + "\t" + chrPos + "\t" +  ref + "\t" +  alt
			
			# Parse info column (ex. DP=##;TI=####;RankSumPos=####;...)	
			infoCol = lineSplit[34].replace('"', '')
			if "SB" in infoCol:
				interestedInfoFields.append("SB")

			infoDict = parseAndFilterInfoValues(infoCol, interestedInfoFields)

			# Parse genotype column (ex. GT:GQ:AD...)
			genotypeHeaderCol = lineSplit[35].replace('"', '')
			genotypeValCol = lineSplit[36].replace('"', '')

			if "SB" in genotypeHeaderCol:
				interestedGenoFields.append("SB")

			genoDict = parseAndFilterGenoValues(genotypeHeaderCol, genotypeValCol, interestedGenoFields)
			
			# Merge all Dictionaries together to return all fields, info fields, and genotype fields from AnnoVar
			avAllFieldDict = dict(list(avFieldDict.items()) + list(infoDict.items()) + list(genoDict.items()))

			if avAllFieldDict["gene"] in interestedGeneList:
				avDict[variantPos] = avAllFieldDict

	infile.close()
	
	return avDict		

def parseAndFilterInfoValues(infoCol, interestedInfoFields):
	
	infoDict = dict.fromkeys(interestedInfoFields, "NA") # Dictionary 
	
	infoArray= infoCol.split(';')

	for idx in enumerate(infoArray):
		num = idx[0]
		value = idx[1]
	
		if "=" in value:
			infoHeader = value.split('=')[0]
			infoVal = value.split('=')[1]
			
			if infoHeader in infoDict:
				infoDict[infoHeader] = infoVal
		
		"""
			"QD" in infoHeader: # quality score divided by total depth. PASS: QD > 2.0
			"BaseQRankSum" in infoHeader:
			"DP" in infoHeader: # Total depth of reads. PASS: DP > 20
			"HaplotypeScore"
			"ReadPosRankSum" in infoHeader: # Chance of it being FP cuz it's end of reads. ReadPosRankSum > -8
			"MQ" in infoHeader: # Mapping Quality. MQ > 30
			"VQSLOD" in infoHeader: # log odd ratios of it being a true variants. 
			"SB" in infoHeader: #Strand bias: higher = higher bias (Flag if value is >=0?)
			"FS" in infoHeader: #Fisher Strand: Indication of strand bias. FS < 60
		"""

	return infoDict

def parseAndFilterGenoValues(genoHeaderCol, genoValCol, interestedGenoFields):
	genoDict = dict.fromkeys(interestedGenoFields, "NA")
	genoHeader = genoHeaderCol.split(':')
	genoVal = genoValCol.split(':')

	if len(genoHeader) == len(genoVal):
		for idx in enumerate(genoHeader):
			num = idx[0]
			header = idx[1]
			if header in genoDict:
				genoDict[header] = genoVal[num]

			"""
			"GT" # 0/0  - Homo reference 0/1 - Hetero 1/1 Homo alternate
			"AD(Ref,Alt)" = genoArray[idx] # Depth per allele - REF depth / ALT depth
			"DP" # Total depth of reads - Flag if DP < 20
			"GQ" # Genotype quality. Phred-scale confidence that GT is true. Flag if < 30
			"SB" # Strand Bias. Flag as strand bias if value is >= 0
			"""
	return genoDict


