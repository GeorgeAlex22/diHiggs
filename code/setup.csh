#!/usr/bin/csh

# source appropriate ROOT and GCC versions
source /cvmfs/sft.cern.ch/lcg/external/gcc/4.9.1/x86_64-slc6-gcc48-opt/setup.csh
source /cvmfs/sft.cern.ch/lcg/releases/ROOT/6.06.06-5e975/x86_64-slc6-gcc49-opt/bin/thisroot.csh

# set Delphes links
setenv PYTHONPATH $HOME/Delphes-3.4.0/python:$PYTHONPATH
setenv LD_LIBRARY_PATH $HOME/Delphes-3.4.0:$LD_LIBRARY_PATH

# set LHAPDF links
setenv PATH $HOME/local/bin:$PATH
setenv LD_LIBRARY_PATH $HOME/local/lib:$LD_LIBRARY_PATH
setenv PYTHONPATH $HOME/local/python3.9/site-packages:$PYTHONPATH
