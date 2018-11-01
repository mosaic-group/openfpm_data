#! /bin/bash

# Make a directory in /tmp/openfpm_data

echo "Build on: $2 with $3"

mkdir /tmp/openfpm_data_$3
mv * .[^.]* /tmp/openfpm_data_$3
mv /tmp/openfpm_data_$3 openfpm_data

mkdir openfpm_data/src/config

#git clone https://git.mpi-cbg.de/openfpm/openfpm_devices.git openfpm_devices
git --version
git -c core.sshCommand="ssh -vvv" clone git@git.mpi-cbg.de:openfpm/openfpm_devices.git openfpm_devices

cd "$1/openfpm_data"

pre_command=""
sh ./autogen.sh
if [ "$2" == "master" ]; then
  options="$options --disable-gpu"
elif [ "$2" == "sbalzarini-mac-15" ]; then
  options="$options --with-libhilbert=$HOME/$4/LIBHILBERT"
fi

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
make

if [ $? -ne 0 ]; then
    curl -X POST --data "payload={\"icon_emoji\": \":jenkins:\", \"username\": \"jenkins\"  , \"attachments\":[{ \"title\":\"Error:\", \"color\": \"#FF0000\", \"text\":\"$2 failed to compile the openfpm_data test $opt_comp \" }] }" https://hooks.slack.com/services/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
    exit 1
fi

