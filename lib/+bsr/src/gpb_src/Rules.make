#############################################################################
# Project directory structure.
#############################################################################

# build scripts directory
build_scripts_dir := scripts/build

# directory layout
bin_dir     := bin
depend_dir  := depend
include_dir := include
src_dir     := src
test_dir    := test
matlab_dir  := matlab
doc_dir     := doc

#############################################################################
# Architecture determination with mex
#############################################################################

# matlab mex file compilation settings - matlab path
MATLABDIR ?= /usr/local/matlab

# osx file compilation settings
MEXEXT := $(shell $(MATLABDIR)/bin/mexext)
ifeq ($(MEXEXT), mexglx)
	MATLAB_ARCH := glnx86
	ADDITIONAL_CXX_LINK :=
	ADDITIONAL_CXX_BUILD :=
endif
ifeq ($(MEXEXT), mexa64)
	MATLAB_ARCH := glnxa64
	ADDITIONAL_CXX_LINK :=
	ADDITIONAL_CXX_BUILD :=
endif
ifeq ($(MEXEXT), mexmaci)
	MATLAB_ARCH := maci
	ADDITIONAL_CXX_LINK := -L/opt/local/lib
	ADDITIONAL_CXX_BUILD := -arch x86 -I/opt/local/include
endif
ifeq ($(MEXEXT), mexmaci64)
	MATLAB_ARCH := maci64
	ADDITIONAL_CXX_LINK := -L/opt/local/lib
	ADDITIONAL_CXX_BUILD := -arch x86_64 -I/opt/local/include
endif

#############################################################################
# C++ compilation settings.
#############################################################################

# compiler
CXX := g++

# compilation settings - libraries to link
CXX_LINK := -ljpeg -lpng $(ADDITIONAL_CXX_LINK)

# compilation settings - warning flags
CXX_WARN_BASIC := -ansi -Wall -Wno-long-long
CXX_WARN_EXTRA := -Wundef -Wpointer-arith -Wold-style-cast \
		-Woverloaded-virtual -Wsign-promo
CXX_WARN  := $(CXX_WARN_BASIC) $(CXX_WARN_EXTRA)

# compilation settings - build flags
#CXX_BUILD := -pthread -fexceptions -fPIC -O3 -rdynamic
#force 64 bit build (e.g. on macos)
CXX_BUILD := -pthread -fexceptions -fPIC -O3 -rdynamic $(ADDITIONAL_CXX_BUILD)

# compilation settings - all flags
CXX_FLAGS := $(CXX_WARN) $(CXX_BUILD)

# compilation settings - linker flags
CXX_LDFLAGS := $(CXX_FLAGS) $(CXX_LINK)

#############################################################################
# Matlab mex file compilation settings (only used if building mex files).
#############################################################################

# matlab mex file compilation settings - include path for mex header files
MEX_INCLUDE_PATH := $(MATLABDIR)/extern/include

# matlab mex file compilation settings - libraries to link
MEX_LINK := $(CXX_LINK) -lmx -lmex -lmat

# matlab mex file compilation settings - warning flags
MEX_WARN := $(CXX_WARN)

# matlab mex file compilation settings - build flags
MEX_BUILD := $(CXX_BUILD) -DMATLAB_MEX_FILE -D_GNU_SOURCE -DNDEBUG

# matlab mex file compilation settings - all flags
MEX_FLAGS := $(MEX_WARN) $(MEX_BUILD)

# matlab mex file compilation settings - linker flags

# osx matlab mex file compilation settings - linker flags
ifeq ($(MEXEXT), mexglx)
	MEX_LDFLAGS := \
		$(MEX_FLAGS) -shared \
		-Wl,--version-script,"$(MATLABDIR)/extern/lib/$(MATLAB_ARCH)/mexFunction.map" \
		-Wl,--rpath-link,"$(MATLABDIR)/bin/$(MATLAB_ARCH)" \
		-L"$(MATLABDIR)/bin/$(MATLAB_ARCH)" $(MEX_LINK)
endif
ifeq ($(MEXEXT), mexa64)
	MEX_LDFLAGS := \
		$(MEX_FLAGS) -shared \
		-Wl,--version-script,"$(MATLABDIR)/extern/lib/$(MATLAB_ARCH)/mexFunction.map" \
		-Wl,--rpath-link,"$(MATLABDIR)/bin/$(MATLAB_ARCH)" \
		-L"$(MATLABDIR)/bin/$(MATLAB_ARCH)" $(MEX_LINK)
endif
ifeq ($(MEXEXT), mexmaci)
	MEX_LDFLAGS := \
	   $(MEX_FLAGS) -bundle \
	   -Wl,-exported_symbols_list,"$(MATLABDIR)/extern/lib/$(MATLAB_ARCH)/mexFunction.map" \
	   -L"$(MATLABDIR)/bin/$(MATLAB_ARCH)" $(MEX_LINK)
endif
ifeq ($(MEXEXT), mexmaci64)
	MEX_LDFLAGS := \
	   $(MEX_FLAGS) -bundle \
	   -Wl,-exported_symbols_list,"$(MATLABDIR)/extern/lib/$(MATLAB_ARCH)/mexFunction.map" \
	   -L"$(MATLABDIR)/bin/$(MATLAB_ARCH)" $(MEX_LINK)
endif

#############################################################################
# Documentation build settings (only used if generating documentation).

# doxygen c++ source code documentation generator
DOXYGEN := doxygen
