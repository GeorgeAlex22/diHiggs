#
#            M A I N
#           =========

info:
	@echo Please type: 'gmake main.x'
	@echo or:          'gmake clean'
	@echo detected host: `hostname`,  os = `uname`

#################################################################
SRC_cc = $(wildcard src/*.cc)
SRC_C += $(wildcard src/*.C)
SRC_cpp += $(wildcard src/*.cpp)


OBJ_cc = $(patsubst src/%.cc, obj/%.o, $(SRC_cc))
OBJ_C = $(patsubst src/%.C, obj/%.o, $(SRC_C))
OBJ_cpp = $(patsubst src/%.cpp, obj/%.o, $(SRC_cpp))
ANAOBJ = $(OBJ_cc) $(OBJ_C) $(OBJ_cpp)


#################################################################


#CFLAGS = -O
CFLAGS = -g -O2 
#CFLAGS = -O -pg

#-I/afs/cern.ch/sw/lcg/external/Boost/1.50.0_python2.6/x86_64-slc6-gcc47-opt/include/boost-1_55

#INC = -Ihh/  -I/afs/cern.ch/user/a/alexandg/local/include
#INC = -Ihh/  -I/afs/cern.ch/user/a/alexandg/local/include/ -I/afs/cern.ch/sw/lcg/external/Boost/1.50.0_python2.6/x86_64-slc6-gcc47-opt/include/boost-1_55
#INC = -Ihh/  -I/afs/cern.ch/user/a/alexandg/local/include/ -I/afs/cern.ch/sw/lcg/external/Boost/1.55.0_python2.7/x86_64-slc6-gcc47-opt/include/boost-1_55
INC = -Ihh/  -I/afs/cern.ch/user/a/alexandg/local/include/ -I/cvmfs/sft.cern.ch/lcg/external/Boost/1.53.0_python2.7/x86_64-slc6-gcc47-opt/include/boost-1_53/

LDFLAGS = -L/afs/cern.ch/user/a/alexandg/local/lib -lLHAPDF 

PWD	:=	$(shell pwd)

#INC = -Ihh/  -I${PWD}/local/include

#LDFLAGS = -L${PWD}/local/lib -lLHAPDF


SYSTEM:= $(shell uname)

ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs)
ROOTGLIBS    := $(shell root-config --glibs)

INC += $(ROOTCFLAGS)
LDFLAGS += $(ROOTGLIBS) $(ROOTLIBS)

##################################################################
ifeq ($(SYSTEM),SunOS)
  FCOMP = f77
  FFLAGS +=
  CPPCOMP = CC
  CFLAGS += -DOS_SOLARIS_2_5 -DOS_NEW_CHECK -DOS_STL_ASSERT -DOS_NO_WSTRING -DOS_NO_ALLOCATORS
  INC    += -I/afs/cern.ch/sw/lhcxx/specific/@sys/ObjectSpace/pro \
            -I/afs/cern.ch/sw/lhcxx/specific/@sys/ObjectSpace/pro/ospace/stl \
            -I/afs/cern.ch/sw/lhcxx/specific/@sys/ObjectSpace/pro/ospace/std \
            -I/afs/cern.ch/sw/lhcxx/specific/@sys/CLHEP/pro/include
  F2CLIB =  -L/usr/lang/SC3.0/lib -lM77 -lF77 -lsunmath
  LINKER = CC -t #-pg
endif
##################################################################
ifeq ($(SYSTEM),Linux)
  FCOMP = g77
  FFLAGS +=
##  ifeq ($(findstring clued0,$(HOST)), clued0)
#    CPPCOMP = g++ # dont forget to "setup root v3_02_07c -q GCC_2_95_2"

#    CPPCOMP = g++ -std=c++11

    CPPCOMP = gcc -Wno-deprecated

    #CPPCOMP = gcc -Wno-deprecated -pthread -m32
    #CFLAGS += --exceptions --thread_safe
    #LINKER = KCC --exceptions --thread_safe

#    LINKER = gcc -lstdc++
#    LINKER = gcc -lstdc++    -I/afs/cern.ch/sw/lcg/external/Boost/1.50.0_python2.6/x86_64-slc6-gcc47-opt/include/boost-1_55
#  LINKER = gcc -lstdc++    -I/afs/cern.ch/sw/lcg/external/Boost/1.55.0_python2.7/x86_64-slc6-gcc47-opt/include/boost-1_55
     LINKER = gcc -lstdc++    -I/cvmfs/sft.cern.ch/lcg/external/Boost/1.53.0_python2.7/x86_64-slc6-gcc47-opt/include/boost-1_53
#    F2CLIB = -L/fnal/ups/gcc/v3_1/lib/
#    F2CLIB = -L/D0/ups/gcc/v3_1/lib/
##  else
##    CPPCOMP = g++
##    CFLAGS += -Df2cFortran -fno-second-underscore \
##              -fno-automatic -fdollar-ok -fno-backslash -finit-local-zero \
##            -fugly-logint -ftypeless-boz -Wall
##    LINKER = g77
##  endif
  F2CLIB = -lg2c   # -lm -lstdc++
endif
##################################################################
ifeq ($(SYSTEM),IRIX64)
  FCOMP = f77
  FFLAGS += -mips3 -n32
  CPPCOMP = KCC
  CFLAGS += -n32 --exceptions --thread_safe
  LINKER = KCC -n32 --exceptions --thread_safe # -v
  F2CLIB += -lm -L/usr/lib32/mips3 -lftn
endif
##################################################################


obj/%.o :: src/%.cc hh/*.h*
	$(CPPCOMP) -c $(CFLAGS) $(INC) $< -o $@

obj/%.o :: src/%.C hh/*.h*
	$(CPPCOMP) -c $(CFLAGS) $(INC) $< -o $@

obj/%.o :: src/%.cpp hh/*.h*
	$(CPPCOMP) -c $(CFLAGS) $(INC) $< -o $@

##############################################################


main.x:$(ANAOBJ)
	$(LINKER) $(LDFLAGS) $(ANAOBJ) \
	-o main.x

dihiggsToyMC.x:$(ANAOBJ)
	$(LINKER) $(LDFLAGS) $(ANAOBJ) \
	-o dihiggsToyMC.x

clean:
	rm -rf obj/*.o
	rm -f *.x




