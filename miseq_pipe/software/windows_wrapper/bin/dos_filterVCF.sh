#!/bin/bash
echo
echo "Starting Cygwin filterVCF script..."
echo

read -r -p "Please drag run folder into terminal window:" response

echo $response

read -r -p "Please drag Manifest (.bed) file into terminal window:" genelist

responsePath="$(/usr/bin/cygpath -u $response)"
genelistPath="$(/usr/bin/cygpath -u $genelist)"

mkdir -p /home/$USER/logs/filterVCF

runName=$(basename $responsePath)
/home/$USER/miseq_pipe/filterVCF.sh $responsePath $genelistPath | tee/home/$USER/logs/filterVCF/$runName.filter.log


