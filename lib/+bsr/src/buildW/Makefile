#flex -P__cfl_ -s -o$@ $<

MATLABDIR ?= /usr/local/matlab
MEXEXT = $(shell $(MATLABDIR)/bin/mexext)

ifeq ($(MEXEXT), mexmaci)
	CXX_FLAGS = -I/opt/local/include -arch x86
	LFLAGS = -L/opt/local/lib
endif
ifeq ($(MEXEXT), mexmaci64)
	CXX_FLAGS = -I/opt/local/include -arch x86_64
	LFLAGS = -L/opt/local/lib
endif

CXX = g++ -g -O3 -Wall
IFLAGS = -Iutil -Igroup -Ipartition 
LFLAGS = -lm -ljpeg $(LFLAGS)
F90FLAGS = -lgfortran -llapack -lblas

#######################################################################################333

UTILOBJS = util/configFileLexer.o util/configure.o util/exception.o util/random.o util/regex.o \
						util/string.o util/timer.o util/util.o util/kmeans.o util/image.o util/message.o util/segmentation.o

OBJS = ic.o affinity.o smatrix.o 

#######################################################################################333

.SUFFIXES : .o .cc 

.cc.o:
	${CXX} -fPIC $(CXX_FLAGS) -c ${IFLAGS} -o $@ -x c++ $<

#######################################################################################333

all: buildW.$(MEXEXT)

#######################################################################################333

buildW.$(MEXEXT): util/libutil.a

#matlab -nodisplay -nojvm -r "mex buildW.cpp -Iutil smatrix.cc ic.cc affinity.cc util/libutil.a; exit"

#######################################################################################

util/libutil.a: ${UTILOBJS} 
	rm -f util/libutil.a && ar cr util/libutil.a $(UTILOBJS) && ranlib util/libutil.a

#######################################################################################333

clean :
	rm -f util/*.a
	rm -f *.o util/*.o
	rm -f *~ util/*~ 
	rm -f *.mex*

