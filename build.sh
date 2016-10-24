#! /bin/bash

source $HOME/.bashrc

echo "$PATH"

# Make a directory in /tmp/openfpm_data

mkdir /tmp/openfpm_io
mv * .[^.]* /tmp/openfpm_io
mv /tmp/openfpm_io openfpm_io

mkdir openfpm_io/src/config

git clone git@git.mpi-cbg.de:openfpm/openfpm_devices.git openfpm_devices
git clone git@git.mpi-cbg.de:openfpm/openfpm_data.git openfpm_data
git clone git@git.mpi-cbg.de:openfpm/openfpm_pdata.git openfpm_pdata
git clone git@git.mpi-cbg.de:openfpm/openfpm_vcluster.git openfpm_vcluster

cd "$1/openfpm_io"

echo "Compiling on $2"

sh ./autogen.sh
if [ "$2" == "master" ]
then
 sh ./configure CXX=mpic++ --with-hdf5=$HOME/HDF5/bin/h5pcc --disable-gpu
elif [ "$2" == "gin" ]
then
 module load gcc/4.8.2
 module load boost/1.54.0
 sh ./configure CXX=mpic++ --with-boost=/sw/apps/boost/1.54.0/ --with-hdf5=$HOME/HDF5/bin/h5pcc
else
 sh ./configure CXX=mpic++ --with-hdf5=$HOME/HDF5/bin/h5pcc
fi
make

if [ $? -ne 0 ]; then
   curl -X POST --data "payload={\"icon_emoji\": \":jenkins:\", \"username\": \"jenkins\"  , \"attachments\":[{ \"title\":\"Error:\", \"color\": \"#FF0000\", \"text\":\"$2 failed to compile the openfpm_io test \" }] }" https://hooks.slack.com/services/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
   exit 1 ; 
fi

