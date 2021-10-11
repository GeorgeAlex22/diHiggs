#!/usr/bin/csh

cd

source /cvmfs/sft.cern.ch/lcg/external/gcc/4.9.1/x86_64-slc6-gcc48-opt/setup.csh
source /cvmfs/sft.cern.ch/lcg/releases/ROOT/6.06.06-5e975/x86_64-slc6-gcc49-opt/bin/thisroot.csh

wget http://cp3.irmp.ucl.ac.be/downloads/Delphes-3.4.0.tar.gz

tar -xvf Delphes-3.4.0.tar.gz
rm Delphes-3.4.0.tar.gz

cd Delphes-3.4.0

make -j 4

# configure env variables
setenv PYTHONPATH $HOME/Delphes-3.4.0/python:$PYTHONPATH
setenv LD_LIBRARY_PATH $HOME/Delphes-3.4.0:$LD_LIBRARYPATH

