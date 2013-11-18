#!/bin/bash
echo
echo "Cygwin BASH Script dos2cywin.sh"
echo

read -r -p "Please drag run folder into terminal window:" response

mkdir -p /home/$USER/miseq_pipe/logs/createReport

/home/$USER/miseq_pipe/createReport.cygwin.sh $(/usr/bin/cygpath -au $response) | tee /home/$USER/miseq_pipe/logs/createReport/$(basename $response).log
sleep 2
