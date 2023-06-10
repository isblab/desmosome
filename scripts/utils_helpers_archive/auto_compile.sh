#!/bin/bash
module load mpi/openmpi-x86_64
if [ -d build ]; then
    rm -r build
    echo CLEANED BUILD DIRECTORY
fi
if [ -d imp ]; then
    yes | rm -r imp
    echo CLEANED OLD IMP DIRECTORY
fi
cp -R ../imp_backup/imp ./
cd imp
./setup_git.py
echo COPIED THE IMP BACKUP FOLDER
cd modules
if [ $# -eq 2 ]; then
    cp -R ../../$1 ./  # copy the custom module folder
    echo COPIED THE CUSTOM MODULE
fi
git clone https://github.com/salilab/imp-sampcon.git sampcon
echo LOADED SAMPCON. PAUSING FOR 60s
sleep 60s
cd ..
cd ..
mkdir build
cd build
cmake ../imp -DCMAKE_CXX_COMPILER=mpic++ -DCMAKE_C_COMPILER=mpicc -DCMAKE_BUILD_TYPE=Release -G Ninja -DCGAL_DIR=/usr/share/CGAL/cmake
echo FINISHED CMAKE. RUNNING THE COMPILATION
if [ $# -eq 2 ]; then
    echo USING $2 CORES
    ninja -j $2 1> compile.log 2> compile_error.log
elif [ $# -eq 1 ]; then
    echo USING $1 CORES
    ninja -j $1 1> compile.log 2> compile_error.log
else
    echo USING 1 CORE AS DEFAULT
    ninja -j 1 1> compile.log 2> compile_error.log
fi
