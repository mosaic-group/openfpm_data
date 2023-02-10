#! /bin/bash

# check if the directory $1/VCDEVEL exist

<<<<<<< HEAD
if [ -d "$1/VCDEVEL" && -f "$1/VCDEVEL/include/Vc/Vc" ]; then
=======
if [ -d "$1/VCDEVEL" -a -f "$1/VCDEVEL/include/Vc/Vc" ]; then
>>>>>>> origin/master
  echo "VCDEVEL already installed"
  exit 0
fi

wget http://ppmcore.mpi-cbg.de/upload/Vc-1.4.1.tar.gz
#rm -rf Vc
tar -xf Vc-1.4.1.tar.gz
cd Vc-1.4.1
mkdir build
cd build
<<<<<<< HEAD
cmake -DCMAKE_INSTALL_PREFIX:PATH=$1/VCDEVEL ..
=======
cmake -DCMAKE_INSTALL_PREFIX:PATH=$1/VCDEVEL -DCMAKE_C_COMPILER=$3 -DCMAKE_CXX_COMPILER=$4 ..
>>>>>>> origin/master
make
make install

