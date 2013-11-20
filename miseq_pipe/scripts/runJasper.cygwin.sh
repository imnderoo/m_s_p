#!/bin/bash

if [ $# -ne 2 ]
then
	echo "Usage: $(basename $0) in_folder sample_id"
	echo " "
	exit 1
fi

inFolder=$1
sampleId=$2

outFolder=$1/pdf

mkdir -p $outFolder

# ====================================
# Run JasperStarter to generate report
# ====================================

# Variable for JASPER_PATH is defined in the parent script createReport.sh
cd $JASPER_PATH

unix2dos $inFolder/${sampleId}_runqc_sorted.csv &>/dev/null
unix2dos $inFolder/${sampleId}_varqc_sorted.csv &>/dev/null
unix2dos $inFolder/${sampleId}_varsummary_sorted.csv &>/dev/null
unix2dos $inFolder/${sampleId}_exonsummary_sorted.csv &>/dev/null

./jasperstarter.exe pr $(cygpath -aw $JASPER_TEMPLATE/RunQC.jrxml) -w -f pdf -o $(cygpath -aw $outFolder) -t csv --data-file $(cygpath -aw $inFolder/${sampleId}_runqc_sorted.csv) --csv-first-row

./jasperstarter.exe pr $(cygpath -aw $JASPER_TEMPLATE/VarQC.jrxml) -w -f pdf -o $(cygpath -aw $outFolder) -t csv --data-file $(cygpath -aw $inFolder/${sampleId}_varqc_sorted.csv) --csv-first-row

./jasperstarter.exe pr $(cygpath -aw $JASPER_TEMPLATE/VarSummary.jrxml) -w -f pdf -o $(cygpath -aw $outFolder) -t csv --data-file $(cygpath -aw $inFolder/${sampleId}_varsummary_sorted.csv) --csv-first-row

./jasperstarter.exe pr $(cygpath -aw $JASPER_TEMPLATE/ExonSummary.jrxml) -w -f pdf -o $(cygpath -aw $outFolder) -t csv --data-file $(cygpath -aw $inFolder/${sampleId}_exonsummary_sorted.csv) --csv-first-row

# ==========================
# Attach header to boxplots
# ==========================

counter=1

for bxp in $inFolder/coverage_plots/*.png
do
        # Update Image Path to attach boxplot to PDF
        
	bxp=$(cygpath -aw $bxp)
	#escapeBXP=$(echo $bxp | sed -e 's/\//\\\//g')
	escapeBXP=$(echo $bxp | sed -e 's/\\/\\\\\\\\/g')

        sed "s/IMAGE_PATH_HERE/$escapeBXP/g" $JASPER_TEMPLATE/BoxplotsMaster.jrxml > $JASPER_TEMPLATE/Boxplots.jrxml
        
	./jasperstarter pr $(cygpath -aw $JASPER_TEMPLATE/Boxplots.jrxml) -w -f pdf -o $(cygpath -aw $outFolder) -t csv --data-file $(cygpath -aw $inFolder/${sampleId}_runqc_sorted.csv) --csv-first-row

	./jasperstarter pr $(cygpath -aw $JASPER_TEMPLATE/Boxplots.jrxml) -w -f pdf -o $(cygpath -aw $outFolder) -t csv --data-file $(cygpath -aw $inFolder/${sampleId}_runqc_sorted.csv) --csv-first-row

        mv $outFolder/Boxplots.pdf $outFolder/Boxplot$counter.pdf
        counter=$((counter+1))
done



