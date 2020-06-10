#! /bin/bash

hostname=$(hostname)
branch=$3

# Make a directory in /tmp/openfpm_data

cd "openfpm_io"

echo "CHECKING MACHINE"
if [ x"$hostname" == x"cifarm-centos-node.mpi-cbg.de"  ]; then
	echo "CENTOS"
        export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$HOME/openfpm_dependencies/openfpm_io/$branch/HDF5/lib"
	export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$HOME/openfpm_dependencies/openfpm_io/$branch/BOOST/lib"
fi

if [ x"$hostname" == x"cifarm-ubuntu-node"  ]; then
        export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$HOME/openfpm_dependencies/openfpm_io/$branch/HDF5/lib"
	export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$HOME/openfpm_dependencies/openfpm_io/$branch/BOOST/lib"
fi

if [ x"$hostname" == x"cifarm-mac-node.mpi-cbg.de"  ]; then
	echo "MACOS X"
        export DYLD_LIBRARY_PATH="$DYLD_LIBRARY_PATH:$HOME/openfpm_dependencies/openfpm_io/$branch/HDF5/lib"
	export DYLD_LIBRARY_PATH="$DYLD_LIBRARY_PATH:$HOME/openfpm_dependencies/openfpm_io/$branch/BOOST/lib"
fi

pwd

./build/src/io


