#! /bin/bash

# Make a directory in /tmp/OpenFPM_data

export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/usr/local/cuda-6.5/lib64"
export PATH="$PATH:/usr/local/cuda-6.5/bin"

mkdir /tmp/openfpm_data
mv * .[^.]* /tmp/openfpm_data
mv /tmp/openfpm_data OpenFPM_data

git clone git@ppmcore.mpi-cbg.de:incardon/openfpm_devices.git OpenFPM_devices

cd $"OpenFPM_data"

sh ./autogen.sh
sh ./configure
make


