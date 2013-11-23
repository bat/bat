# Root variables
ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := -lMinuit $(shell root-config --libs)

#ROOTLIBS     += -lRooFitCore -lRooFit -lRooStats -lFoam -lMathMore

# compiler and flags
CXX          = g++
CXXFLAGS     = -g -O2 -Wall -fPIC -Wno-deprecated
LD           = @LD@
LDFLAGS      = -g -O2

# standard commands
RM           = rm -f
MV           = mv
ECHO         = echo

# add ROOT flags
CXXFLAGS    += $(ROOTCFLAGS)

# ----------------------------------------------------------------------
# The following definitions depend on the setup of the system where
# the project is being compiled. If BAT is installed in the standard
# system search path or the installation directory is defined in the
# BATINSTALLDIR environmental variable then the lines below are correct
# and the compilation will work
CXXFLAGS    += -I. -I./include -I$(BATINSTALLDIR)/include
LIBS        += -L$(BATINSTALLDIR)/lib -lBATmodels -lBATmtf -lBAT $(ROOTLIBS)

#LIBS        +=  -lcuba

# List of all classes (models) used in the program
# Add classes to the end. Backslash indicates continuation
# on the next line
CXXSRCS      = \
	BCPostProcessor.cxx \
	BCPPMarginalize.cxx \
	BCPPDiagnostics.cxx

# ----------------------------------------------------------------------
# don't change lines below unless you know what you're doing
#

CXXOBJS      = $(patsubst %.cxx,%.o,$(CXXSRCS))
EXEOBJS      =
MYPROGS     = \
	runPlotter \
	runMarginalize

GARBAGE      = $(CXXOBJS) $(EXEOBJS) *.o *~ link.d $(MYPROGS)


# targets
all : plotter marginalize

link.d : $(patsubst %.cxx,%.h,$(CXXSRCS))
	$(CXX) -MM $(CXXFLAGS) $(CXXSRCS) > link.d;

-include link.d

%.o : %.cxx
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean :
	$(RM) $(GARBAGE)

plotter : runPlotter.cxx $(CXXOBJS)
	$(CXX) $(CXXFLAGS) -c $<
	$(CXX) $(LDFLAGS) runPlotter.o $(CXXOBJS) $(LIBS) -o runPlotter

marginalize : runMarginalize.cxx $(CXXOBJS)
	$(CXX) $(CXXFLAGS) -c $<
	$(CXX) $(LDFLAGS) runMarginalize.o $(CXXOBJS) $(LIBS) -o runMarginalize

print :
	echo compiler  : $(CXX)
	echo c++ srcs  : $(CXXSRCS)
	echo c++ objs  : $(CXXOBJS)
	echo c++ flags : $(CXXFLAGS)
	echo libs      : $(LIBS)
	echo so flags  : $(SOFLAGS)

	echo rootlibs  : $(ROOTLIBS)

