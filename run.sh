#! /bin/bash

hostname=$(hostname)
branch=$3

# Make a directory in /tmp/openfpm_data

cd "$1/openfpm_io"

if [ "$2" == "gin" ]
then
 module load gcc/4.9.2
 module load boost/1.54.0
fi

if [ x"$hostname" == x"cifarm-centos-node.mpi-cbg.de"  ]; then
        export LD_LiBRARY_PATH="$LD_LIBRARY_PATH:$HOME/openfpm_dependencies/openfpm_io/$branch/HDF5/lib"
	ls $HOME/openfpm_dependencies/openfpm_io/$branch/HDF5/lib
fi

if [ x"$hostname" == x"cifarm-ubuntu-node"  ]; then
        export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$HOME/openfpm_dependencies/openfpm_io/$branch/HDF5/lib"
	ls $HOME/openfpm_dependencies/openfpm_io/$branch/HDF5/lib
fi

if [ x"$hostname" == x"cifarm-mac-node.mpi-cbg.de"  ]; then
        export DYLD_LIBRARY_PATH="$DYLD_LIBRARY_PATH:$HOME/openfpm_dependencies/openfpm_io/$branch/HDF5/lib"
	ls $HOME/openfpm_dependencies/openfpm_io/$branch/HDF5/lib
fi

./build/src/io
if [ $? -ne 0 ]; then
   curl -X POST --data "payload={\"icon_emoji\": \":jenkins:\", \"username\": \"jenkins\"  , \"attachments\":[{ \"title\":\"Error:\", \"color\": \"#FF0000\", \"text\":\"$2 failed to complete the openfpm_io test \" }] }" https://hooks.slack.com/services/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
   exit 1 ; 
fi


