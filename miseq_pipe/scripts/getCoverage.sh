#!/bin/bash

if [ $# -ne 4 ]
then
     	echo "----------------------------------------------------------------------------------------------------------------"
     	echo "This script creates coverage plots and stats for the input BAM file using the regions found in BED file."
	echo "----------------------------------------------------------------------------------------------------------------"
	echo "Usage: nohup $(basename $0) input_folder mapq_threshold  gene_list region_of_interest.bed &"
	echo "     nohup makes sure script continues to even when terminal closes."
	echo "     & makes sure script runs in background)"
	echo " "
	echo "Example: nohup ./$(basename $0) ./NGS_data/Run_3/ 20 ./genelists/george.breastcancer.180413.genelist ./genelists/george.breastcancer.180413.bed &"
	echo " "
	exit 1
fi

input_folder=$1

mkdir -p $input_folder/coverage/

echo "$(basename $0) started for $1"
echo " "

qScore=$2
genelist=$(readlink -f $3)
bedfile=$(readlink -f $4)

for bam in $(find $1 -name "*q$qScore.bam")
do
	# echo "$bam $genelist $bedfile"
	coverage.sge $bam $genelist $bedfile 15
done

echo " "
echo "$(basename $0) finished for $1"
