#! /bin/bash

# Make a directory in /tmp/OpenFPM_data

mkdir /tmp/openfpm_data
mv * .[^.]* /tmp/openfpm_data
mv /tmp/openfpm_data OpenFPM_data

git clone git@ppmcore.mpi-cbg.de:incardon/openfpm_devices.git OpenFPM_devices

cd "$1/OpenFPM_data"

sh ./autogen.sh
sh ./configure
make


