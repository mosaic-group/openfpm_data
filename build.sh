#! /bin/bash

# Make a directory in /tmp/openfpm_data

echo "Build on: $2 with $4 branch $5"

# Check if libHilbert is installed

if [ ! -d $HOME/openfpm_dependencies/openfpm_data/LIBHILBERT ]; then
	./install_LIBHILBERT.sh $HOME/openfpm_dependencies/openfpm_data/ 4
fi

echo "AAAAAAAAAAAAAAAAAAA $2"
if [ "$2" == x"mac" ]; then
	echo "REMOVING"
	rm -rf $HOME/openfpm_dependencies/openfpm_data/BOOST
fi

if [ ! -d $HOME/openfpm_dependencies/openfpm_data/BOOST ]; then
	if [ "$2" == x"mac" ]; then
        	./install_BOOST.sh $HOME/openfpm_dependencies/openfpm_data/ 4 darwin
	else
		./install_BOOST.sh $HOME/openfpm_dependencies/openfpm_data 4 gcc
	fi
fi

ls $HOME/openfpm_dependencies/openfpm_data/BOOST/lib

mkdir /tmp/openfpm_data_$3
mv * .[^.]* /tmp/openfpm_data_$3
mv /tmp/openfpm_data_$3 openfpm_data

mkdir openfpm_data/src/config

git clone https://git.mpi-cbg.de/openfpm/openfpm_devices.git openfpm_devices
cd openfpm_devices
git checkout GPU_test
cd ..

cd "$1/openfpm_data"

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
    curl -X POST --data "payload={\"icon_emoji\": \":jenkins:\", \"username\": \"jenkins\"  , \"attachments\":[{ \"title\":\"Error:\", \"color\": \"#FF0000\", \"text\":\"$2 failed to comfigure openfpm_data test $opt_comp \" }] }" https://hooks.slack.com/services/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
    exit 1
fi
make VERBOSE=1

if [ $? -ne 0 ]; then
    curl -X POST --data "payload={\"icon_emoji\": \":jenkins:\", \"username\": \"jenkins\"  , \"attachments\":[{ \"title\":\"Error:\", \"color\": \"#FF0000\", \"text\":\"$2 failed to compile the openfpm_data test $opt_comp \" }] }" https://hooks.slack.com/services/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
    exit 1
fi

