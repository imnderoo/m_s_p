#!/bin/bash

if [ $# -ne 2 ]
then
     	echo "----------------------------------------------------------------------------------------------------------------"
     	echo "This script filters out aligned reads in the BAM file that have a mapping quality score less than specified threshold"
	echo "----------------------------------------------------------------------------------------------------------------"
	echo "Usage: nohup $(basename $0) input_folder mapq_threshold &"
	echo "     nohup makes sure script continues to even when terminal closes."
	echo "     & makes sure script runs in background)"
	echo " "
	echo "Example: nohup ./$(basename $0) ./NGS_data/Test_Run/ 30 &"
	echo " "
	exit 1
fi

input_folder=$1

echo "$(basename $0) started for $1"
echo " "

qScore=$2

# =====================
# Remove PCR Duplicates
# =====================

for bam in $(find $1 -name "*.bam")
do
	echo "Removing PCR Duplicates... for $bam"
	if [[ $bam =~ .*dedup*.* ]]
	then 
		echo "$bam already deduped"
	
	elif [[ $bam =~ .*q\d*.* ]]
	then
		echo "$bam already filtered using qScore"
	else
		input_bam=$(basename $bam)
		output_bam=$(echo $input_bam | sed 's/.bam/_dedup.bam/')
	
		samtools rmdup $bam $input_folder/$output_bam &>/dev/null &

		pid_last_rmdup=$!
		sleep 8
	fi
done

wait $pid_last_rmdup

# ========================
# Filter using MAPQ Scores
# ========================

for bam in $(find $1 -name "*_dedup.bam")
do
	echo "Filtering by MAPQ scores for $bam ..."
	if [[ $bam =~ .*q\d*.* ]]
	then
		echo "$bam already filtered using qScore"
	else
		input_bam=$(basename $bam)
		output_bam=$(echo $input_bam | sed 's/.bam/_q'$qScore'.bam/') 
                
		samtools view -bh -F d -q$qScore $bam > $input_folder/$output_bam &
		pidf_last_filter=$!
		sleep 15
	fi
done

wait $pid_last_filter

echo " "
echo "$(basename $0) finished for $1"
