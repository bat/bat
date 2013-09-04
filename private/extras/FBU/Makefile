# Root variables
ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs) -lMinuit
ROOTGLIBS    := $(shell root-config --glibs)

# directories
INCLUDEDIR = include
SRCDIR     = src
OBJDIR     = obj

# compiler and flags
CXX          = g++
CXXFLAGS     = -g -Wall -fPIC -Wno-deprecated -O2
LD           = g++
LDFLAGS      = -g -O2 -Wall
SOFLAGS      = -shared

# directories
INCLUDEDIR = include
SRCDIR     = src
OBJDIR     = obj

# standard commands
RM           = rm -f
MV           = mv
ECHO         = echo
CINT         = rootcint

# add ROOT flags
CXXFLAGS    += $(ROOTCFLAGS) 
CXXFLAGS    += -I. -I./include -I$(BATINSTALLDIR)/include
LIBS        += $(ROOTLIBS)  -L$(BATINSTALLDIR)/lib -lBAT -lBATmodels
GLIBS       += $(ROOTGLIBS) -L$(BATINSTALLDIR)/lib -lBAT -lBATmodels

# List of all classes (models) used in the program
# Add classes to the end. Backslash indicates continuation
# on the next line
CXSRCS      = \
        BCFBU.cxx \
        BCFBUBackground.cxx \
	BCFBUSystematic.cxx \
	BCFBUNormSystematic.cxx \
	BCFBUBkgSystematic.cxx

# ----------------------------------------------------------------------
# don't change lines below unless you know what you're doing
#


HEADERSRCS = $(patsubst %.cxx, $(INCLUDEDIR)/%.h, $(CXSRCS))
CXXSRCS = $(patsubst %.cxx, $(SRCDIR)/%.cxx, $(CXSRCS))
CXXOBJS = $(patsubst %.cxx, $(OBJDIR)/%.o, $(CXSRCS))
EXEOBJS =
MYPROGS = runBCMultiTemplateFitter
LIBSO   = libFBU.so
GARBAGE = $(CXXOBJS) $(EXEOBJS) *.o *~ link.d $(MYPROGS) $(LIBSO)

# targets
all : library

link.d : $(patsubst %.cxx,%.h,$(CXXSRCS))
	$(CXX) -MM $(CXXFLAGS) $(CXXSRCS) > link.d;

-include link.d

%.o : %.cxx
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean :
	$(RM) $(GARBAGE)

library : $(CXXOBJS)
	$(CXX) $(CXXFLAGS) $(LIBS) $(SOFLAGS) $^ -o $(LIBSO)

# rule for single class output files
$(OBJDIR)/%.o : $(SRCDIR)/%.cxx $(INCLUDEDIR)/%.h $(PREPROCESSOR_H)
	@if [ ! -e $(OBJDIR) ]; then $(MKDIR) $(OBJDIR); fi
	$(CXX) $(CXXFLAGS) -c $< -o $@

tarball : 
	tar -cf MTF.tar \
	Makefile \
	include/*.h \
	src/*.cxx \
	examples/test_SB/CreateHistograms.C \
	examples/test_SB/*.cxx \
	examples/test_SB/Makefile \
	examples/test_SB/setup.sh \
	examples/test_SB_ensemble1/CreateHistograms.C \
	examples/test_SB_ensemble1/*.cxx \
	examples/test_SB_ensemble1/Makefile \
	examples/test_SB_ensemble1/setup.sh \
	examples/test_SB_ensemble2/CreateHistograms.C \
	examples/test_SB_ensemble2/*.cxx \
	examples/test_SB_ensemble2/Makefile \
	examples/test_SB_ensemble2/setup.sh \
	examples/test_SB_calibration/CreateHistograms.C \
	examples/test_SB_calibration/*.cxx \
	examples/test_SB_calibration/Makefile \
	examples/test_SB_calibration/setup.sh \
	examples/test_SB_channels/CreateHistograms.C \
	examples/test_SB_channels/*.cxx \
	examples/test_SB_channels/Makefile \
	examples/test_SB_channels/setup.sh \
	examples/test_SB_function/CreateHistograms.C \
	examples/test_SB_function/*.cxx \
	examples/test_SB_function/Makefile \
	examples/test_SB_function/setup.sh \
	examples/test_SB_templatefunctions/CreateHistograms.C \
	examples/test_SB_templatefunctions/*.cxx \
	examples/test_SB_templatefunctions/Makefile \
	examples/test_SB_templatefunctions/setup.sh \
	examples/test_whel/CreateHistograms.C \
	examples/gerda/*.cxx \
	examples/gerda/Makefile \
	examples/gerda/setup.sh \
	examples/test_whel/CreateHistograms.C \
	examples/test_whel/*.cxx \
	examples/test_whel/Makefile \
	examples/test_whel/setup.sh \
	doc/manual_MTF.pdf
	gzip MTF.tar

print :
	echo compiler  : $(CXX)
	echo c++ srcs  : $(CXXSRCS)
	echo c++ objs  : $(CXXOBJS)
	echo c++ flags : $(CXXFLAGS)
	echo libs      : $(LIBS)
	echo so flags  : $(SOFLAGS)

	echo rootlibs  : $(ROOTLIBS)
	echo rootglibs : $(ROOTGLIBS)

