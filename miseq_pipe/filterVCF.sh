#!/bin/bash

if [ $# -ne 2 ]
then
        echo "-----------------------------------------------------------------------------------------------------------------"
        echo "This script goes through all the .bam and .vcf files for each sample and generates a report found in NGS_analysis"
        echo "-----------------------------------------------------------------------------------------------------------------"
        echo "Usage: nohup $(basename $0) input_folder genelists/manifest.bed &"
        echo "     nohup makes sure script continues to even when terminal closes."
        echo "     & makes sure script runs in background)"
        echo " "
        echo "Example: nohup ./$(basename $0) ./NGS_data/Run_3/ &"
        echo " "
        exit 1
fi

# TODO: Before continuing, check that the required files and folders exist.
# TODO: Namely "Summary.xml" and Interop folder"

# =========================
# Variables
# =========================
runFolder=$1
genelist=$2
#genelist=$HOME/miseq_pipe/genelists/george.breastcancer.120613.bed

vcftoolsPath=$HOME/miseq_pipe/software/vcftools_0.1.7
export PATH=$PATH:$vcftoolsPath

rename .vcf .vcf.original *
find $runFolder -name "*.vcf.original" -print0 | xargs -0 -I file vcftools --vcf file --recode --bed $genelist --out file
rename .vcf.original.recode.vcf .vcf *
