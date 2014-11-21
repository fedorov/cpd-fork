CC=g++
LD=g++

DEBUG = 1

BUILDDIR = ./build
SRCDIR = ./
TET_DIR = ../../ThirdParty/tetgen1.4.3/source
INCLUDEDIR = $(TET_DIR)
BINDIR = ./

CPPFLAGS=-I$(INCLUDEDIR) --std=c++11
ifdef DEBUG
CPPFLAGS:=$(CPPFLAGS) -g -O0 -DPTR_EXCEEDS_LONG -DTETLIBRARY
else
CPPFLAGS:=$(CPPFLAGS) -O3
endif
LDFLAGS = 

# Library extension
ifndef LIB_EXT
LIB_EXT:=dll
endif

# Matlab folders
ifndef MATLAB_ROOT
MATLAB_ROOT:="C:/Program Files/MATLAB/R2013b"
# MATLAB_ROOT:="/Applications/MATLAB_R2013b.app"
endif
ifndef MATLAB_BINDIR
MATLAB_BINDIR:="$(MATLAB_ROOT)/bin/win64"
# MATLAB_BINDIR:="$(MATLAB_ROOT)/bin/maci64"
endif
ifndef MEX_EXT
MEX_EXT:=mexw64
# MEX_EXT:=mexmaci64
endif
ifndef MATLAB_INCLUDEDIR
MATLAB_INCLUDEDIR:="$(MATLAB_ROOT)/extern/include"
endif

MEX_CPPFLAGS:=$(CPPFLAGS) -DMX_COMPAT_32 -DMATLAB_MEX_FILE -I"$(MATLAB_INCLUDEDIR)" -Wall
MEX_LDFLAGS:= $(LDFLAGS) -shared -L"$(MATLAB_BINDIR)" -lstdc++ -lmex -lmx -lmat

MEX_SOURCES = tetgen_mex.cpp
MEX_OBJECTS = $(BUILDDIR)/tetgen_mex.mexo
MEX = tetgen_mex.$(MEX_EXT)

SOURCES := $(TET_DIR)/predicates.cxx $(TET_DIR)/tetgen.cxx 
OBJECTS := $(BUILDDIR)/predicates.o $(BUILDDIR)/tetgen.o

MKDIR_CMD=mkdir -p $(@D)

default: all

cleanup:
	rm -rf $(OBJECTS) $(MEX_OBJECTS)

clean: cleanup
	rm -rf $(MEX) 

# do not delete intermediates
.SECONDARY:

vars:
	@echo "MEX_SOURCES: $(MEX_SOURCES)"
	@echo "MEX_OBJECTS: $(MEX_OBJECTS)"
	@echo "MEX_CPPFLAGS: $(MEX_CPPFLAGS)"
	@echo "MEX_LDFLAGS: $(MEX_LDFLAGS)"
	@echo "MEX: $(MEX)"
	@echo "SOURCES: $(SOURCES)"
	@echo "OBJECTS: $(OBJECTS)"
	

all: mex

mex: $(MEX)

$(BUILDDIR)/%.o: $(TET_DIR)/%.cxx
	@$(MKDIR_CMD)
	@echo Compiling $@...
	$(CC) $(CPPFLAGS) -o $@ -c $<

$(BUILDDIR)/%.mexo: $(SRCDIR)/%.cpp
	@$(MKDIR_CMD)
	@echo Compiling $@...
	@$(CC) $(MEX_CPPFLAGS) -o $@ -c $<

tetgen_mex.$(MEX_EXT): $(MEX_OBJECTS) $(OBJECTS)
	@$(MKDIR_CMD)
	@echo Assembling $@ ...
	$(LD) -o $@ $(MEX_LDFLAGS) $(OBJECTS) $<
	@echo Complete!
	
