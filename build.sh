#! /bin/bash

# Make a directory in /tmp/openfpm_data

workspace=$1
hostname=$(hostname)
type_compile=$3
branch=$4

echo "Build on: $hostname with $type_compile branch: $branch"

# Check if libHilbert is installed

if [ ! -d $HOME/openfpm_dependencies/openfpm_data/LIBHILBERT ]; then
	./install_LIBHILBERT.sh $HOME/openfpm_dependencies/openfpm_data/ 4
fi


if [ ! -d $HOME/openfpm_dependencies/openfpm_data/BOOST ]; then
	if [ x"$hostname" == x"cifarm-mac-node" ]; then
		echo "Compiling for OSX"
        	./install_BOOST.sh $HOME/openfpm_dependencies/openfpm_data/ 4 darwin
	else
		echo "Compiling for Linux"
		./install_BOOST.sh $HOME/openfpm_dependencies/openfpm_data 4 gcc
	fi
fi

mkdir /tmp/openfpm_data_$type_compile
mv * .[^.]* /tmp/openfpm_data_$type_compile
mv /tmp/openfpm_data_$type_compile openfpm_data

mkdir openfpm_data/src/config

git clone https://git.mpi-cbg.de/openfpm/openfpm_devices.git openfpm_devices
cd openfpm_devices
git checkout GPU_test
cd ..

cd "$workspace/openfpm_data"

pre_command=""
sh ./autogen.sh
options="$options --disable-gpu "
options="$options --with-boost=$HOME/openfpm_dependencies/openfpm_data/BOOST  --with-libhilbert=$HOME/openfpm_dependencies/openfpm_data/LIBHILBERT"

if [ x"$3" == x"SE"  ]; then
  options="$options --enable-se-class1 --enable-se-class2 --enable-se-class3 --with-action-on-error=throw --enable-test-coverage"
  opt_comp="for security enhancement"
fi

if [ x"$3" == x"VALGRIND" ]; then
  pre_command="valgrind --leak-check=full"
  options="$options --disable-gpu --enable-test-coverage"
  opt_comp="for valgrind test"
fi

sh ./configure $options
if [ $? -ne 0 ]; then
    curl -X POST --data "payload={\"icon_emoji\": \":jenkins:\", \"username\": \"jenkins\"  , \"attachments\":[{ \"title\":\"Error:\", \"color\": \"#FF0000\", \"text\":\"$hostname failed to comfigure openfpm_data test $opt_comp \" }] }" https://hooks.slack.com/services/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
    exit 1
fi
make

if [ $? -ne 0 ]; then
    curl -X POST --data "payload={\"icon_emoji\": \":jenkins:\", \"username\": \"jenkins\"  , \"attachments\":[{ \"title\":\"Error:\", \"color\": \"#FF0000\", \"text\":\"$hostname failed to compile the openfpm_data test $opt_comp \" }] }" https://hooks.slack.com/services/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
    exit 1
fi

ls build/src/mem_map


