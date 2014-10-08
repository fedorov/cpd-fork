###############################################################################
#                                                                             #
# makefile for TetGen                                                         #
#                                                                             #
# Type "make" to compile TetGen into an executable program (tetgen).          #
# Type "make tetlib" to compile TetGen into a library (libtet.a).             #
# Type "make distclean" to delete all object (*.o) files.                     #
#                                                                             #
###############################################################################

# CXX should be set to the name of your favorite C++ compiler.
# ===========================================================

CXX = g++
#CXX = icpc
#CXX = CC

# CXXFLAGS is the level of optimiztion, default is -O. One should try
# -O2, -O3 ... to find the best optimization level.
# ===================================================================

CXXFLAGS = -g

# PREDCXXFLAGS is for compiling J. Shewchuk's predicates. It should
# always be equal to -O0 (no optimization). Otherwise, TetGen may not
# work properly.

PREDCXXFLAGS = -O0

# SWITCHES is a list of switches to compile TetGen.
# =================================================
#
# By default, TetGen uses double precision floating point numbers.  If you
#   prefer single precision, use the -DSINGLE switch. 
#
# The source code of TetGen includes a lot of assertions, which are mainly
#   used for catching bugs at that places.  These assertions somewhat slow
#   down the speed of TetGen.  They can be skipped by define the -DNDEBUG
#   switch.

SWITCHES = -Wall -DSELF_CHECK

# SWITCHES = -Wall -Wabi -Wctor-dtor-privacy \
#            -Woverloaded-virtual -Wno-pmf-conversions -Wsign-promo \
#            -Wsynth  -Wchar-subscripts -Wconversion -Wsign-compare \
#            -Wcomment  -Wimplicit -Wmissing-braces -Wparentheses \
#            -Wreturn-type -Wswitch -Wswitch-default \
#            -Wswitch-enum -Wtrigraphs -W -DSELF_CHECK

# RM should be set to the name of your favorite rm (file deletion program).
RM = /bin/rm -f
MKDIR_P = mkdir -p

# CYGWIN PTR BUG

MACHINE = $(shell uname -m) # not sure how portable this is ...
ifeq ($(findstring 64,$(MACHINE)),64)
	ARCH = 64
	CXXFLAGS += -DPTR_EXCEEDS_LONG -DAMD_64
	PREDCXXFLAGS += -DPTR_EXCEEDS_LONG -DAMD_64
else
	ARCH = 32
endif


# INPUT/OUTPUT DIRECTORIES

# change these to set the proper directories where each files shoould be
SRCDIR   = .
OBJDIR   = ../build
BINDIR   = ../build/Release
LIBDIR   = ../build/Release

SOURCES  := $(wildcard $(SRCDIR)/*.cxx)
INCLUDES := $(wildcard $(SRCDIR)/*.h)
OBJECTS  := $(SOURCES:$(SRCDIR)/%.cxx=$(OBJDIR)/%.o)

# The action starts here.


default: tetgen

all: tetgen tetlib

tetgen:	$(SRCDIR)/tetgen.cxx $(OBJDIR)/predicates.o
	@$(MKDIR_P) $(BINDIR)
	$(CXX) $(CXXFLAGS) $(SWITCHES) -o $(BINDIR)/tetgen $(SRCDIR)/tetgen.cxx $(OBJDIR)/predicates.o -lm

tetlib: $(SRCDIR)/tetgen.cxx $(OBJDIR)/predicates.o
	@$(MKDIR_P) $(LIBDIR)
	$(CXX) $(CXXFLAGS) $(SWITCHES) -DTETLIBRARY -o $(OBJDIR)/tetgen.o -c $<
	ar r $(LIBDIR)/libtet.a $(OBJDIR)/tetgen.o $(OBJDIR)/predicates.o

$(OBJDIR)/predicates.o: $(SRCDIR)/predicates.cxx
	@$(MKDIR_P) $(OBJDIR)
	$(CXX) $(PREDCXXFLAGS) -c $< -o $@

.PHONEY: clean
clean:
	@$(RM) $(OBJECTS)
	@echo "Cleanup complete!"

.PHONEY: remove
remove: clean
	@$(RM) $(BINDIR)/tetgen $(LIBDIR)/libtet.a
	@echo "Executable removed!"




