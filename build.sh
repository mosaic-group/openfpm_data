#! /bin/bash


workspace=$1
hostname=$(hostname)
branch=$3

echo "Directory: $workspace"
echo "Machine: $hostname"
echo "Branch name: $branch"

echo "Branch: $3"

# Make a directory in /tmp/openfpm_data

mkdir /tmp/openfpm_io
mv * .[^.]* /tmp/openfpm_io
mv /tmp/openfpm_io openfpm_io

git clone git@git.mpi-cbg.de:openfpm/openfpm_devices.git openfpm_devices
cd openfpm_devices
git checkout master
cd ..
git clone git@git.mpi-cbg.de:openfpm/openfpm_data.git openfpm_data
cd openfpm_data
git checkout master
cd ..
git clone git@git.mpi-cbg.de:openfpm/openfpm_pdata.git openfpm_pdata
cd openfpm_pdata
git checkout master
cd ..
git clone git@git.mpi-cbg.de:openfpm/openfpm_vcluster.git openfpm_vcluster
cd openfpm_vcluster
git checkout master
cd ..

cd "$1/openfpm_io"

rm -rf $HOME/openfpm_dependencies/openfpm_io/$branch/MPI
rm -rf $HOME/openfpm_dependencies/openfpm_io/$branch/HDF5
rm -rf $HOME/openfpm_dependencies/openfpm_io/$branch/BOOST

if [ x"$hostname" == x"cifarm-centos-node.mpi-cbg.de"  ]; then
	source /opt/rh/devtoolset-7/enable
        ./install_MPI_mpich.sh $HOME/openfpm_dependencies/openfpm_io/$branch/ 4
	export PATH="/opt/bin:$HOME/openfpm_dependencies/openfpm_io/$branch/MPI/bin/:$PATH"
        ./install_BOOST.sh $HOME/openfpm_dependencies/openfpm_io/$branch/ 4
	./install_HDF5.sh $HOME/openfpm_dependencies/openfpm_io/$branch/ 4
fi

if [ x"$hostname" == x"cifarm-ubuntu-node"  ]; then
        ./install_MPI_mpich.sh $HOME/openfpm_dependencies/openfpm_io/$branch/ 4
	export PATH="$HOME/openfpm_dependencies/openfpm_io/$branch/MPI/bin/:$PATH"
	./install_BOOST.sh $HOME/openfpm_dependencies/openfpm_io/$branch/ 4
	./install_HDF5.sh $HOME/openfpm_dependencies/openfpm_io/$branch/ 4
fi

if [ x"$hostname" == x"cifarm-mac-node.mpi-cbg.de"  ]; then
        export PATH="/usr/local/bin:$PATH"
        ./install_MPI_mpich.sh $HOME/openfpm_dependencies/openfpm_io/$branch/ 4
	export PATH="$HOME/openfpm_dependencies/openfpm_io/$branch/MPI/bin/:$PATH"
	./install_BOOST.sh $HOME/openfpm_dependencies/openfpm_io/$branch/ 4
	./install_HDF5.sh $HOME/openfpm_dependencies/openfpm_io/$branch/ 4
fi

# Go in the right branch

echo "Compiling on $2"

sh ./autogen.sh
sh ./configure CXX=mpic++ --with-hdf5=$HOME/openfpm_dependencies/openfpm_io/$branch/HDF5 --with-boost=$HOME/openfpm_dependencies/openfpm_io/$branch/BOOST --with-pdata=../../openfpm_pdata/ --enable-cuda-on-cpu

make VERBOSE=1 -j 4


