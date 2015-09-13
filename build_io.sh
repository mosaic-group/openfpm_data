#! /bin/bash

# Make a directory in /tmp/openfpm_data

mkdir /tmp/openfpm_io
mv * .[^.]* /tmp/openfpm_io
mv /tmp/openfpm_io openfpm_io

mkdir openfpm_io/src/config

git clone git@ppmcore.mpi-cbg.de:incardon/openfpm_devices.git openfpm_devices
git clone git@ppmcore.mpi-cbg.de:incardon/openfpm_data.git openfpm_data
git clone git@ppmcore.mpi-cbg.de:incardon/openfpm_pdata.git openfpm_pdata
cd "$1/openfpm_data"
git checkout develop
cd ..

cd "$1/openfpm_IO"

echo "Compiling on $2"

sh ./autogen.sh
if [ "$2" == "master" ]
then
 sh ./configure --disable-gpu
elif [ "$2" == "gin" ]
then
 module load gcc/4.8.2
 module load boost/1.54.0
 sh ./configure --with-boost=/sw/apps/boost/1.54.0/
else
 sh ./configure
fi
make

./src/io

