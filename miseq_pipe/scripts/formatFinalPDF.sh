#!/bin/bash

if [ $# -ne 1 ]
then
	echo "Usage: $(basename $0) <sample_folder_path>"
	echo " "
	exit 1
fi

inFolder=$(readlink -f $1)
sampleid=$(basename $inFolder)

# ============================
# Merge Jasper Reports into One
# ============================

boxplots=""

# Find all boxplots path
for bxp in $inFolder/pdf/Boxplot*.pdf
do
	boxplots="$boxplots $bxp"
done

# Merge using PDFTK

mergeOutFile=$inFolder/pdf/${sampleid}_merge.pdf

pdftk \
$inFolder/pdf/VarSummary.pdf \
$inFolder/pdf/ExonSummary.pdf \
$inFolder/pdf/VarQC.pdf \
$inFolder/pdf/RunQC.pdf \
$boxplots \
cat output \
$mergeOutFile

# ============================
# Add Page Numbers to the PDF Report
# ============================

escapeReport=$(echo $mergeOutFile | sed -e 's/\//\\\//g')

texOutFile=$inFolder/pdf/${sampleid}_final.tex
pdfOutFile=$inFolder/pdf/${sampleid}_final.pdf

sed "s/PDF_PATH_HERE/$escapeReport/g" $LATEX_TEMPLATE > $texOutFile

# Switch to the pdf directory so that the output tex and pdf files are in the pdf folder
cd $inFolder/pdf/

pdflatex -draftmode $texOutFile &>/dev/null 
pdflatex $texOutFile &>/dev/null

mv $pdfOutFile $inFolder
