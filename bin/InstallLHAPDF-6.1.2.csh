 ## Set up LCG build tools (optional) and install directory

source /afs/cern.ch/sw/lcg/contrib/gcc/4.7/x86_64-slc6-gcc47-opt/setup.csh

mkdir local


wget http://www.hepforge.org/archive/lhapdf/LHAPDF-6.1.2.tar.gz
tar -xvf LHAPDF-6.1.2.tar.gz

cd LHAPDF-6.1.2

#./configure --prefix=$HOME/local

#./configure --prefix=$HOME/local --with-boost=/afs/cern.ch/sw/lcg/external/Boost/1.50.0_python2.6/x86_64-slc6-gcc47-opt/include/boost-1_55

./configure --prefix=$HOME/local --with-boost=/afs/cern.ch/sw/lcg/external/Boost/1.55.0_python2.7/x86_64-slc6-gcc47-opt/include/boost-1_55/

make -j2 && make install

cd $HOME/local/share/LHAPDF

# Currently, these two wget commands do not seem to fetch a file that can be extracted. 
# I downloaded the tarballs locally and then moved the to lxplus with scp
# (e.g. scp -r Downloads/CT10[nlo].tar.gz <username>@lxplus.cern.ch:~local/share/LHAPDF/CT10[nlo].tar.gz)

wget http://www.hepforge.org/archive/lhapdf/pdfsets/6.1/CT10nlo.tar.gz 
tar -xvf CT10nlo.tar.gz
rm -r CT10nlo.tar.gz

wget http://www.hepforge.org/archive/lhapdf/pdfsets/6.1/CT10.tar.gz 
tar -xvf CT10.tar.gz
rm -r CT10.tar.gz


cd ..


setenv PATH $PWD/local/bin:$PATH
setenv LD_LIBRARY_PATH $PWD/local/lib:$LD_LIBRARY_PATH
setenv PYTHONPATH $PWD/local/lib64/python2.6/site-packages:$PYTHONPATH

lhapdf-config --help

lhapdf list