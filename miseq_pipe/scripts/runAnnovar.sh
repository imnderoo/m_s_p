#!/bin/bash

if [ $# -ne 3 ]
then
     	echo "----------------------------------------------------------------------------------------------------------------"
	echo "This script runs annovar on all the VCF files found in input_folder. The output is found in input_folder/annovar"
     	echo "----------------------------------------------------------------------------------------------------------------"
	echo "Usage: nohup ./$(basename $0) input_folder bed genelist variantdb &"
	echo "     nohup makes sure script continues to even when terminal closes."
	echo "     & makes sure script runs in background)"
	echo " "
	echo "Example: nohup ./$(basename $0) ./NGS_data/Test_Run ./genelists/george.breastcancer.180413.genelist ./variantdb/george.breastcancer.format.db &"
	echo " "
	exit 1
fi

input_folder=$1
genelist=$2
variantsDB=$3

name=$(echo $input_folder | sed 's/\/$//' | awk -F/ '{print $NF}')
echo $name

mkdir -p $input_folder/annovar/intermediate/
mkdir -p $input_folder/annovar/input/

echo "runAnnovar.sh started for $input_folder"

#===========================
# Run AnnoVar
#===========================

for vcfs in $(ls $input_folder/*.vcf)
do
	if [[ $vcfs =~ *filtered* ]]
	then
		echo "$vcfs already filtered"
	else
		#filtered_vcfs=$(echo $vcfs | sed 's/.vcf/_filtered.vcf/')
		#echo "Filtering VCF using BED File $vcfs"
		#java -Xmx2g -jar $GATK_TOOL -R $REF -T SelectVariants --variant $vcfs -o $filtered_vcfs -L $roi_bed
		
		echo "Running AnnoVAR for $vcfs"
		annovar.sge $vcfs vcf
	fi
done

#=============================
# Organize output into folders
#=============================

mv $input_folder/annovar/*annovar*.* $input_folder/annovar/intermediate/
mv $input_folder/annovar/*annovar* $input_folder/annovar/input/
mv $input_folder/annovar/intermediate/*multianno.txt $input_folder/annovar/

echo "runAnnovar.sh finished for $input_folder"
