#! /bin/bash

# Make a directory in /tmp/OpenFPM_data

mkdir /tmp/openfpm_io
mv * .[^.]* /tmp/openfpm_io
mv /tmp/openfpm_io OpenFPM_IO

mkdir OpenFPM_IO/src/config

git clone git@ppmcore.mpi-cbg.de:incardon/openfpm_devices.git OpenFPM_devices
git clone git@ppmcore.mpi-cbg.de:incardon/openfpm_data.git OpenFPM_data

cd "$1/OpenFPM_IO"

sh ./autogen.sh
if [ "$2" == "master" ]
then
 sh ./configure --disable-gpu
else
 sh ./configure
fi
make

./src/io

