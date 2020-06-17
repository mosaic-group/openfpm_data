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
git checkout sparse_cl
cd ..
git clone git@git.mpi-cbg.de:openfpm/openfpm_data.git openfpm_data
cd openfpm_data
git checkout sparse_cl_broken
cd ..
git clone git@git.mpi-cbg.de:openfpm/openfpm_pdata.git openfpm_pdata
cd openfpm_pdata
git checkout sparse_cl
cd ..
git clone git@git.mpi-cbg.de:openfpm/openfpm_vcluster.git openfpm_vcluster
cd openfpm_vcluster
git checkout sparse_cl
cd ..

cd "$1/openfpm_io"


if [ x"$hostname" == x"cifarm-centos-node.mpi-cbg.de"  ]; then
        ./install_MPI_mpich.sh $HOME/openfpm_dependencies/openfpm_io/$branch/ 4
	export PATH="$HOME/openfpm_dependencies/openfpm_io/$branch/MPI/bin/:$PATH"
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
if [ "$2" == "master" ]
then
 sh ./configure CXX=mpic++ --with-hdf5=$HOME/$3/HDF5/bin/h5pcc --disable-gpu
elif [ "$2" == "gin" ]
then
 module load gcc/4.8.2
 module load boost/1.54.0
 sh ./configure CXX=mpic++ --with-boost=/sw/apps/boost/1.54.0/ --with-hdf5=$HOME/$3/HDF5/bin/h5pcc
else
 sh ./configure CXX=mpic++ --with-hdf5=$HOME/openfpm_dependencies/openfpm_io/$branch/HDF5 --with-boost=$HOME/openfpm_dependencies/openfpm_io/$branch/BOOST --with-pdata=../../openfpm_pdata/
fi

make VERBOSE=1 -j 4

if [ $? -ne 0 ]; then
   curl -X POST --data "payload={\"icon_emoji\": \":jenkins:\", \"username\": \"jenkins\"  , \"attachments\":[{ \"title\":\"Error:\", \"color\": \"#FF0000\", \"text\":\"$2 failed to compile the openfpm_io test \" }] }" https://hooks.slack.com/services/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
   exit 1 ; 
fi

