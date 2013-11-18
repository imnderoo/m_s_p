#!/bin/bash

if [ $# -ne 3 ]
then
	echo " "
	echo "Usage: $(basename $0) <input_file> <output_file> <gene_list>"	
	echo " "
	echo "The Input file should be the variants MS Access database table exported as a XLS file."
	echo "The exported XLS file should have had the columns PubMedString, DNA Reference and DNA Interpretation removed"
	echo "Finally, the XLS file should be saved as a tab-delimited .txt file."
	exit 1
fi

input=$1
output=$2
gene=$3

perl -pi -e "s/\r/\n/g" $input


perl /home/george/scripts/format_variant_db.pl $input $output $gene
