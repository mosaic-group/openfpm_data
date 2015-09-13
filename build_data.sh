#! /bin/bash

# Make a directory in /tmp/OpenFPM_data

mkdir /tmp/openfpm_data
mv * .[^.]* /tmp/openfpm_data
mv /tmp/openfpm_data openfpm_data

mkdir OpenFPM_data/src/config

git clone git@ppmcore.mpi-cbg.de:incardon/openfpm_devices.git openfpm_devices

cd "$1/OpenFPM_data"

sh ./autogen.sh
if [ "$2" == "master" ]
then
 sh ./configure --disable-gpu
else
 sh ./configure
fi
make

./src/mem_map

