###################################################################
# This Makefile was created using the ./CreateProject.sh script
# for project StackTool
# ./CreateProject.sh is part of Bayesian Analysis Toolkit (BAT).
# BAT can be downloaded from http://www.mppmu.mpg.de/bat
###################################################################
#
# Run 'make' to compile the program and 'make clean' to remove
# all compiled parts and 'clean' the directory.
#
# You might need to adjust the CFLAGS, LIBS, and GLIBS based on
# the BAT installation on your system. Consult the gmake manual
# for details.
#
###################################################################

# Root variables
ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs) -lMinuit
ROOTGLIBS    := $(shell root-config --glibs)

# compiler and flags
CXX          = g++
CXXFLAGS     = -g -Wall -fPIC -Wno-deprecated -O2
LD           = g++
LDFLAGS      = -g -O2
SOFLAGS      = -shared

# standard commands
RM           = rm -f
MV           = mv
ECHO         = echo
CINT         = rootcint

# add ROOT flags
CXXFLAGS    += $(ROOTCFLAGS)

# ----------------------------------------------------------------------
# The following definitions depend on the setup of the system where
# the project is being compiled. If BAT is installed in the standard
# system search path then the lines below are correct and the compilation
# will work
CXXFLAGS    += -I. -I./include
LIBS        += $(ROOTLIBS)  -L../SummaryTool -lSummaryTool
GLIBS       += $(ROOTGLIBS) -L../SummaryTool -lSummaryTool

# In case you see following errors during the compilation,
#
#   undefined reference to 'Divonne'
#   undefined reference to 'Suave'
#   undefined reference to 'Cuhre'
#   undefined reference to 'Vegas'
#
# your version of BAT was installed with the Cuba support and you need
# to adjust the Makefile by uncommenting the following lines. You might
# also need to add the path to libcuba.a to the lines below
# as -L/path/to/cuba/lib
#
# LIBS        += -lcuba
# GLIBS       += -lcuba

# If BAT was installed in a non-standard path (e.g. user home
# directory) then one can specify the location of the header
# files here (uncomment following lines and adjust the path in
# BATINSTALLDIR):
#
# BATINSTALLDIR = '/path/to/bat/installation/directory'
 CXXFLAGS    += -I$(BATINSTALLDIR)/include
 LIBS        += -L$(BATINSTALLDIR)/lib -lBAT 
 GLIBS       += -L$(BATINSTALLDIR)/lib -lBAT

# List of all classes (models) used in the program
# Add classes to the end. Baskslash indicates continuation
# on the next line
CXXSRCS      = \
        StackModel.cxx \
        StackModelManager.cxx \
				EnsembleTestTool.cxx

# ----------------------------------------------------------------------
# don't change lines below unless you know what you're doing
#

CXXOBJS      = $(patsubst %.cxx,%.o,$(CXXSRCS))

GARBAGE      = $(CXXOBJS) $(EXEOBJS) *.o *~ link.d $(MYPROGS) libStackTool.so

# targets
all : libstacktool

link.d : $(patsubst %.cxx,%.h,$(CXXSRCS))
	$(CXX) -MM $(CXXFLAGS) $(CXXSRCS) > link.d;

-include link.d

%.o : %.cxx
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean :
	$(RM) $(GARBAGE)

hypothesistest : runHypothesisTest.cxx $(CXXOBJS)
	$(CXX) $(CXXFLAGS) -c $<
	$(CXX) $(LDFLAGS) $(LIBS) runHypothesisTest.o $(CXXOBJS) -o runHypothesisTest

templatefit : runTemplateFit.cxx $(CXXOBJS)
	$(CXX) $(CXXFLAGS) -c $<
	$(CXX) $(LDFLAGS) $(LIBS) runTemplateFit.o $(CXXOBJS) -o runTemplateFit

ensembletest : runEnsembleTest.cxx $(CXXOBJS)
	$(CXX) $(CXXFLAGS) -c $<
	$(CXX) $(LDFLAGS) $(LIBS) runEnsembleTest.o $(CXXOBJS) -o runEnsembleTest

gerda : runGerda.cxx $(CXXOBJS)
	$(CXX) $(CXXFLAGS) -c $<
	$(CXX) $(LDFLAGS) $(LIBS) runGerda.o $(CXXOBJS) -o runGerda

templategerda : runTemplateGerda.cxx $(CXXOBJS)
	$(CXX) $(CXXFLAGS) -c $<
	$(CXX) $(LDFLAGS) $(LIBS) runTemplateGerda.o $(CXXOBJS) -o runTemplateGerda

libstacktool: $(CXXOBJS)
	$(CXX) $(SOFLAGS) $(LDFLAGS) $^ -o libStackTool.so

print :
	echo compiler  : $(CXX)
	echo c++ srcs  : $(CXXSRCS)
	echo c++ objs  : $(CXXOBJS)
	echo c++ flags : $(CXXFLAGS)
	echo libs      : $(LIBS)
	echo so flags  : $(SOFLAGS)

	echo rootlibs  : $(ROOTLIBS)
	echo rootglibs : $(ROOTGLIBS)

tarball :
	mkdir StackTool_1.0.0
	mkdir StackTool_1.0.0/Whel
	mkdir StackTool_1.0.0/hypothesis
	cp README StackTool_1.0.0/.
	cp setup.sh StackTool_1.0.0/.
	cp Makefile StackTool_1.0.0/.
	cp StackModel.h StackTool_1.0.0/.
	cp StackModel.cxx StackTool_1.0.0/.
	cp StackModelManager.h StackTool_1.0.0/.
	cp StackModelManager.cxx StackTool_1.0.0/.
	cp EnsembleTestTool.h StackTool_1.0.0/.
	cp EnsembleTestTool.cxx StackTool_1.0.0/.
	cp Whel/runTemplateFit.cxx StackTool_1.0.0/Whel/.
	cp Whel/runEnsembleTest.cxx StackTool_1.0.0/Whel/.
	cp Whel/runTemplateFit.cxx StackTool_1.0.0/Whel/.
	cp Whel/Makefile StackTool_1.0.0/Whel/.
	cp hypothesis/runHypothesisTest.cxx StackTool_1.0.0/hypothesis/.
	cp hypothesis/histograms.root StackTool_1.0.0/hypothesis/.
	cp hypothesis/Makefile StackTool_1.0.0/hypothesis/.
	tar -cf StackTool_1.0.0.tar StackTool_1.0.0
	gzip StackTool_1.0.0.tar
	rm -r StackTool_1.0.0