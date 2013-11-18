#!/bin/bash

if [ $# -ne 4 ]
then
     	echo "----------------------------------------------------------------------------------------------------------------"
	echo "This script runs annovar on all the VCF files found in input_folder. The output is found in input_folder/annovar"
     	echo "----------------------------------------------------------------------------------------------------------------"
	echo "Usage: nohup ./$(basename $0) input_folder bed genelist variantdb &"
	echo "     nohup makes sure script continues to even when terminal closes."
	echo "     & makes sure script runs in background)"
	echo " "
	echo "Example: nohup ./$(basename $0) ./NGS_data/Test_Run ./genelists/george.breastcancer.180413.genelist ./variantdb/george.breastcancer.format.db ./genelists/george.breastcancer010913.exoncoords  &"
	echo " "
	exit 1
fi

input_folder=$1
genelist=$2
variantsDB=$3
geneExons=$4

name=$(echo $input_folder | sed 's/\/$//' | awk -F/ '{print $NF}')

mkdir -p $input_folder/annovar/intermediate/
mkdir -p $input_folder/annovar/input/

#===========================
# Flatten AnnoVar Columns
#===========================

echo "Flattening AnnoVAR"

for avs in $(ls $input_folder/annovar/*multianno.txt)
do
	flatten_annovar.py $avs "DP" "GT,GQ,AD,SB"
done

#===========================
# Filter AnnoVar Fields
#===========================

echo "Filtering AnnoVAR"

for avs in $(ls $input_folder/annovar/*multianno_flatten.txt)
do

filter_annovar.py $avs $genelist $variantsDB $geneExons \
"Chr,Start,End,Ref,Alt,\
Func.refGene,Gene.refGene,ExonicFunc.refGene,AAChange.refGene,\
esp6500si_all,1000g2012apr_all,snp137,\
Zygosity,QUAL,FILTER,DP,AD,GQ,GT,SB"

done
