#! /bin/bash

source $HOME/.bashrc
source $HOME/openfpm_vars_$3

echo "$PATH"

# Make a directory in /tmp/openfpm_data

cd "$1/openfpm_io"

if [ "$2" == "gin" ]
then
 module load gcc/4.8.2
 module load boost/1.54.0
fi

./src/io
if [ $? -ne 0 ]; then
   curl -X POST --data "payload={\"icon_emoji\": \":jenkins:\", \"username\": \"jenkins\"  , \"attachments\":[{ \"title\":\"Error:\", \"color\": \"#FF0000\", \"text\":\"$2 failed to complete the openfpm_io test \" }] }" https://hooks.slack.com/services/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
   exit 1 ; 
fi


