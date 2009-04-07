# Root variables
ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs) -lMinuit
ROOTGLIBS    := $(shell root-config --glibs)

# Programs
CXX          = g++
CXXFLAGS     = -g -Wall -fPIC -Wno-deprecated -O2
LD           = g++
LDFLAGS      = -g -O2
SOFLAGS      = -shared

RM           = rm -f
MV           = mv
ECHO         = echo
CINT         = rootcint

# Assign or Add variables
CXXFLAGS    += $(ROOTCFLAGS)
CXXFLAGS    += -I. -I./include -I$(BATINSTALLDIR)/BAT -I$(ROOTSYS)/roofit/roofitcore/inc -I$(ROOTSYS)/roofit/roofit/inc

LIBS        += $(ROOTLIBS)  -L$(BATINSTALLDIR)/lib -l RooFitCore -l RooFit -lBAT
GLIBS       += $(ROOTGLIBS) -L$(BATINSTALLDIR)/lib -l RooFitCore -l RooFit -lBAT

CXSRCS      =         BCRooInterface.cxx

CXXSRCS      = $(patsubst %.cxx,%.cxx,$(CXSRCS))
CXXOBJS      = $(patsubst %.cxx,%.o,$(CXSRCS))
EXEOBJS      =
GARBAGE      = $(CXXOBJS) $(EXEOBJS) *.o *~ link.d BCRooInterface.exe

all : project

link.d : $(patsubst %.cxx,%.h,$(CXSRCS))
	$(CXX) -MM $(CXXFLAGS) $(CXXSRCS) > link.d;

include link.d

%.o : %.cxx
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean :
	$(RM) $(GARBAGE)

project : runRooInterface.cxx $(CXXOBJS)
	$(CXX) $(CXXFLAGS) -c $<
		$(CXX) $(LDFLAGS) $(LIBS) runRooInterface.o \
		$(CXXOBJS) -o runRooInterface.exe

print :
	echo compiler  : $(CXX)
	echo compiler  : $(CXSRCS)
	echo c++ srcs  : $(CXXSRCS)
	echo c++ objs  : $(CXXOBJS)
	echo c++ flags : $(CXXFLAGS)
	echo libs      : $(LIBS)
	echo so flags  : $(SOFLAGS)

	echo rootlibs  : $(ROOTLIBS)
	echo rootglibs : $(ROOTGLIBS)

