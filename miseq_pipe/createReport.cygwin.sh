#!/bin/bash


function usage
{
        echo "-----------------------------------------------------------------------------------------------------------------"
        echo "This script goes through all the .bam and .vcf files for each sample and generates a report found in NGS_analysis"
        echo "-----------------------------------------------------------------------------------------------------------------"
        echo "Usage: nohup $(basename $0) input_folder [options] &"
        echo "     nohup makes sure script continues to even when terminal closes."
        echo "     & makes sure script runs in background)"
        echo " "
        echo "Example: nohup ./$(basename $0) ./NGS_data/Run_3/ &"
        echo " "
	echo "Options"
	echo "-------"
	echo "--report-only	Skips Coverage and AnnoVar to generate report. Good for testiny."
	exit 1
}


if [ $# -lt 1 ]
then
	usage
fi

# ================
# Parse Parameters
# ================
runFolder=$1
alignFolder="$1/Data/Intensities/BaseCalls/Alignment"
analysisRootFolder="$HOME/NGS_analysis"

while [ "$2" != "" ]
do
	case $2	in
		--report-only )		report_only="true"
	esac
	shift
done

if [ -n "$report_only" ]
then
	echo "Option --report-only Set"
else
	echo "$report_only"
fi

function checkRequiredFiles
{

	if [ ! -d "$alignFolder" ]
	then
		echo "Error: The required folder '$alignFolder' is not found."
		echo "Please make sure that the desired alignment folder is named Alignment"
		exit 1
	fi
	if [ ! -d "$runFolder/InterOp" ]
	then
		echo "Error: The required folder 'InterOp' is not found in $runFolder"
		exit 1
	fi
	if [ ! -e "$runFolder/RunInfo.xml" ]
	then
		echo "Error: The required file 'RunInfo.xml' is not found in $runFolder"
		exit 1
	fi
	if [ ! -e "$alignFolder/Summary.xml" ]
	then
		echo "Error: The required file 'Summary.xml' is not found in $alignFolder"
		exit 1
	fi
	if [ ! -e "$runFolder/SampleSheet.csv" ]
	then
		echo "Error: The required file 'SampleSheet.csv' is not found in $runFolder"
		exit 1
	fi
}

checkRequiredFiles $runFolder $alignFolder

START_TIME=$(date +%y.%m.%d\|%H:%M:%S)

# =====================================
  echo "SETTING UP VARIABLES" 
# ======================================
# Required Tools: samtools, jasperstarter, coverageBed
createReportFolder="$HOME/miseq_pipe"

samtoolsPath="$HOME/miseq_pipe/software/samtools-0.1.19/bin"
bedtoolsPath="$HOME/miseq_pipe/software/bedtools-2.17.0/bin"
pdftexPath="$HOME/miseq_pipe/software/pdftex-1.40.10/build-pdftex/texk/web2c"
export AVLOC="$createReportFolder/software/annovar"

# Required Scripts
scriptPath="$createReportFolder/scripts"
helperScriptsPath="$createReportFolder/scripts/helpers"
pythonPath="$createReportFolder/scripts/python_modules"

# Export to path
export PATH=$PATH:$scriptPath:$helperScriptsPath:$samtoolsPath:$bedtoolsPath:$pdftexPath
export PYTHONPATH=$pythonPath

# Required Report Templates
export LATEX_TEMPLATE="$createReportFolder/report/addPages.cygwin.tex"
export JASPER_TEMPLATE="$createReportFolder/report/template"
export JASPER_PATH="$HOME/miseq_pipe/software/jasperstarter/bin"
export REFGENE="$createReportFolder/genelists/refGene.broad.sorted.txt"

# Reference Files and Constants
qscore=30
genelist="$createReportFolder/genelists/george.breastcancer.180413.genelist"
genebed="$createReportFolder/genelists/george.breastcancer.120613.bed"
geneExonCoords="$createReportFolder/genelists/george.breastcancer.010913.exoncoords"
variantsdb="$createReportFolder/variantdb/george.breastcancer.300813.txt"
annovarOutputSuffix="_annovar_vcf.hg19_multianno"
coverageOutputSuffix="_dedup_q$qscore"

# =========================
# Fix up VariantsDB 
# =========================
sed -i -e 's/\r\n/\n/g' $variantsdb
sed -i -e 's/\r/\n/g' $variantsdb

ngslib=$(basename $runFolder)
echo "NGSLIB: $ngslib"

if [ ! -n "$report_only" ]
then
	# =========================
	  echo "FILTERING BAM FILE AND REMOVING DUPLICATES"
	# =========================

	filterBAM.sh $alignFolder $qscore

	# ==========================
	  echo "CALCULATING COVERAGE"
	# ==========================

	getCoverage.sh $alignFolder $qscore $genelist $genebed &

	pid_coverage=$!


	# =========================
	  echo "RUNNING ANNOVAR"
	# =========================
	
	AV_START_TIME=$(date +%y.%m.%d\|%H:%M:%S)
	avlog=$HOME/miseq_pipe/logs/av
	mkdir -p $avlog

	runAnnovar.sh $alignFolder $genelist $variantsdb &>$avlog/av_$ngslib.log
	
	wait $pid_coverage
	
	AV_END_TIME=$(date +%y.%m.%d\|%H:%M:%S)
	
	echo "annoVar started at: $AV_START_TIME ended at: $AV_END_TIME"

fi

# =========================
# Filter AnnoVar
# =========================
filterAnnovar.sh $alignFolder $genelist $variantsdb $geneExonCoords

# =========================
# Create Output Directories
# =========================
mkdir -p $analysisRootFolder/$ngslib/Reports

for samples in $(find $alignFolder -name "*.vcf")
do
	sampleid=$(basename $samples | rev | cut -d. -f2- | rev)

	outFolder=$analysisRootFolder/$ngslib/$sampleid

	echo "Generating Report for $sampleid"
	
	mkdir -p $outFolder/coverage_plots
	mkdir -p $outFolder/temp
	mkdir -p $outFolder/csv
	mkdir -p $outFolder/annovar


	# Copy AV Raw files
	cp $alignFolder/annovar/$sampleid${annovarOutputSuffix}.txt $outFolder/annovar/${sampleid}_av.txt
	cp $alignFolder/annovar/$sampleid${annovarOutputSuffix}_flatten.txt $outFolder/annovar/${sampleid}_av_flatten.txt

	# Copy AV files to PDF info files
	cp $alignFolder/annovar/$sampleid${annovarOutputSuffix}_flatten_filter.csv $outFolder/${sampleid}_varqc.tmp
	cp $alignFolder/annovar/$sampleid${annovarOutputSuffix}_flatten_filter.csv $outFolder/${sampleid}_varsummary.tmp
	# Copy Coverage files
	cp $alignFolder/coverage/$sampleid${coverageOutputSuffix}_resequence.csv $outFolder/${sampleid}_exonsummary.tmp
	cp $alignFolder/coverage/gene_coverage/boxplot/$sampleid$coverageOutputSuffix*.png $outFolder/coverage_plots
	# Copy RunQC files
	report_qc_miseq.py $runFolder $outFolder/${sampleid}_varqc.tmp $(basename $variantsdb) "./$(basename $0) $runFolder"	

	# Copy BAM files
	cp $alignFolder/$sampleid.bam $outFolder/
	cp $alignFolder/$sampleid.bam.bai $outFolder/

	# Sort the .tmp files using the correct columns
	sort_csv.py $outFolder/${sampleid}_varsummary.tmp ReportCat,Gene.refGene,ACMG
	sort_csv.py $outFolder/${sampleid}_varqc.tmp Gene.refGene
	sort_csv.py $outFolder/${sampleid}_exonsummary.tmp filter,gene

	# Concatenate SampleID,FlowCell,RunDate to the front of the .tmp and save as .csv
	add_run_info_to_csv.py $alignFolder $outFolder/${sampleid}_varqc_sorted.tmp
	add_run_info_to_csv.py $alignFolder $outFolder/${sampleid}_varsummary_sorted.tmp
	add_run_info_to_csv.py $alignFolder $outFolder/${sampleid}_exonsummary_sorted.tmp
	
	echo "Running JasperReports"	

	# Use Jasper Reports to compile the PDF reports
	runJasper.cygwin.sh $outFolder $sampleid

	# Merge and Add Page Numbers to the PDF report
	formatFinalPDF.sh $outFolder

	# Move temp files to temp folder
	mv $outFolder/*.tmp $outFolder/temp
	mv $outFolder/*.csv $outFolder/csv
	echo "Reports for $sampleid are saved in $outFolder"

	cp $outFolder/${sampleid}_final.pdf $analysisRootFolder/$ngslib/Reports
done


echo "Reports completed"

END_TIME=$(date +%y.%m.%d\|%H:%M:%S)

echo "Scripted started at $START_TIME and ended at $END_TIME"
