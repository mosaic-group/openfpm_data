#! /bin/bash

# Make a directory in /tmp/openfpm_data

mkdir /tmp/openfpm_data
mv * .[^.]* /tmp/openfpm_data
mv /tmp/openfpm_data openfpm_data

mkdir openfpm_data/src/config

git clone git@ppmcore.mpi-cbg.de:incardon/openfpm_devices.git openfpm_devices

cd "$1/openfpm_data"

sh ./autogen.sh
if [ "$2" == "master" ]
then
 sh ./configure --disable-gpu
else
 sh ./configure
fi
make

./src/mem_map

