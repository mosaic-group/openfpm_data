#! /bin/bash

# Make a directory in /tmp/openfpm_data

echo "Build on: $2 with $3"

mkdir /tmp/openfpm_data_$3
mv * .[^.]* /tmp/openfpm_data_$3
mv /tmp/openfpm_data_$3 openfpm_data

mkdir openfpm_data/src/config

git clone git@ppmcore.mpi-cbg.de:incardon/openfpm_devices.git openfpm_devices

cd "$1/openfpm_data"

pre_command=""
sh ./autogen.sh
if [ "$2" == "master" ]; then
  options="$options --disable-gpu"
fi

if [ x"$3" == x"SE"  ]; then
  options="$options --enable-se-class1 --enable-se-class2 --enable-se-class3 --with-action-on-error=stop --enable-test-coverage"
fi

if [ x"$3" == x"VALGRIND" ]; then
  precommand="valgrind"
  options="$options --disable-gpu"
fi

sh ./configure $options
if [ $? -ne 0 ]; then
    exit 1
fi
make
if [ $? -ne 0 ]; then
    exit 1
fi

$pre_command ./src/mem_map
if [ $? -ne 0 ]; then
    exit 1
fi

