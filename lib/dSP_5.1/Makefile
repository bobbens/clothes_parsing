#update the path variables according to your needs and choose CFLAGS to compile respective versions of the programs

#MATLABHOME	= /site/opt/matlab/r2012b/x64
#MATLABHOME = /opt/matlab
#MATLABHOME  =/usr/local/MATLAB/R2012a
MATLABHOME = /usr/local/matlab-R2012a
MEXSUFFIX	= mexa64
CC	= g++
MCC	= mpic++
CFLAGS	= -std=c++98 -pedantic -W -Wall -fopenmp -O3 -fPIC
#CFLAGS	+= -DWITH_FLOAT
#CFLAGS	+= -DWITH_PARALLELINF
#CFLAGS	+= -DWITH_NOPARALLEL
#CFLAGS	+= -DWITH_LATENTCLUSTER

#no changes should be required below this line

MEXFLAGS	= -cxx CC='$(CC)' CXX='$(CC)' LD='$(CC)'
MEXFLAGSMPI	= -cxx CC='$(MCC)' CXX='$(MCC)' LD='$(MCC)'
MEX	= $(MATLABHOME)/bin/mex
INCLUDE	= -I$(MATLABHOME)/extern/include

all: prog_dSP prog_mdSP prog_dLSP structuredPrediction.$(MEXSUFFIX) latentStructuredPrediction.$(MEXSUFFIX)

clean:
	rm -f *.a
	rm -f *.o
	rm -f *.$(MEXSUFFIX)
	rm -f prog_dSP
	rm -f prog_mdSP
	rm -f prog_dLSP

libRegionBPP.a: libRegionBPP/libRegionBPP.h libRegionBPP/libRegionBPP.cpp libRegionBPP/DataContainers.h GeneralDefinitions.h
	$(CC) libRegionBPP/libRegionBPP.cpp -c $(CFLAGS)
	ar rcs libRegionBPP.a libRegionBPP.o

libSPRegion.a: libSPRegion/libSPRegion.h libSPRegion/libSPRegion.cpp libSPRegion/Instantiation.cpp libSPRegion/libSPBase.h libSPRegion/libSPRegionMultM.h libSPRegion/libSPRegionMultM.cpp GeneralDefinitions.h
	$(MCC) libSPRegion/libSPRegion.cpp libSPRegion/Instantiation.cpp libSPRegion/libSPRegionMultM.cpp -c $(CFLAGS)
	ar rcs libSPRegion.a libSPRegion.o Instantiation.o libSPRegionMultM.o

libDH.a: libDH/libDHBase.h libDH/libDHStandard.h libDH/libDHStandard.cpp GeneralDefinitions.h
	$(CC) libDH/libDHStandard.cpp -c $(CFLAGS)
	ar rcs libDH.a libDHStandard.o

libDHMultM.a: libDHMultM/libDHBaseMultM.h libDHMultM/libDHSampPMultM.h libDHMultM/libDHSampPMultM.cpp libDHMultM/libDHModlPMultM.h libDHMultM/libDHModlPMultM.cpp GeneralDefinitions.h
	$(MCC) libDHMultM/libDHSampPMultM.cpp libDHMultM/libDHModlPMultM.cpp -c $(CFLAGS)
	ar rcs libDHMultM.a libDHSampPMultM.o libDHModlPMultM.o

libLatentSPRegion.a: libLatentSPRegion/libLatentSPRegion.h libLatentSPRegion/libLatentSPRegion.cpp GeneralDefinitions.h
	$(MCC) libLatentSPRegion/libLatentSPRegion.cpp -c $(CFLAGS)
	ar rcs libLatentSPRegion.a libLatentSPRegion.o

prog_dSP: Test_dSP/Test_dSP.cpp GeneralDefinitions.h libRegionBPP.a libSPRegion.a libDH.a
	$(CC) Test_dSP/Test_dSP.cpp -o prog_dSP $(CFLAGS) -L. -lSPRegion -lDH -lRegionBPP

prog_mdSP: Test_mdSP/Test_mdSP.cpp GeneralDefinitions.h libRegionBPP.a libSPRegion.a libDHMultM.a
	$(MCC) Test_mdSP/Test_mdSP.cpp -o prog_mdSP $(CFLAGS) -L. -lSPRegion -lDHMultM -lRegionBPP

prog_dLSP: Test_dLSP/Test_dLSP.cpp GeneralDefinitions.h libRegionBPP.a libSPRegion.a libDH.a libDHMultM.a libLatentSPRegion.a
	$(MCC) Test_dLSP/Test_dLSP.cpp -o prog_dLSP $(CFLAGS) -L. -lLatentSPRegion -lSPRegion -lDHMultM -lDH -lRegionBPP

structuredPrediction.o: structuredPrediction/structuredPrediction.cpp GeneralDefinitions.h
	$(CC) $(CFLAGS) -static $(INCLUDE) -c structuredPrediction/structuredPrediction.cpp

latentStructuredPrediction.o: structuredPrediction/structuredPrediction.cpp GeneralDefinitions.h
	$(MCC) $(CFLAGS) -DWITH_LATENT -static $(INCLUDE) -c structuredPrediction/structuredPrediction.cpp -o latentStructuredPrediction.o

structuredPrediction.$(MEXSUFFIX): structuredPrediction.o libRegionBPP.a libSPRegion.a libDH.a
	$(MEX) -v $(MEXFLAGS) -output structuredPrediction structuredPrediction.o -L. -lSPRegion -lDH -lRegionBPP -lgomp

latentStructuredPrediction.$(MEXSUFFIX): latentStructuredPrediction.o libRegionBPP.a libSPRegion.a libDH.a libLatentSPRegion.a libDHMultM.a
	$(MEX) -v $(MEXFLAGSMPI) -output latentStructuredPrediction latentStructuredPrediction.o -L. -lLatentSPRegion -lSPRegion -lDH -lDHMultM -lRegionBPP -lgomp
