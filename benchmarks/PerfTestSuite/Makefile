# Root variables
ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs) -lMinuit 
ROOTGLIBS    := $(shell root-config --glibs)

# Programs
CXX          = g++
CXXFLAGS     = -g -Wall -fPIC -Wno-deprecated -O0
LD           = g++
LDFLAGS      = -g -O0
SOFLAGS      = -shared

AR           = ar
ARFLAGS      = rc
RANLIB       = ranlib

RM           = rm -f
MV           = mv
ECHO         = echo
CINT         = rootcint

# Assign or Add variables
CXXFLAGS    += $(ROOTCFLAGS)
CXXFLAGS    += -I./include -I./src -I./
LIBS        += $(ROOTLIBS) -L. -L$(BATINSTALLDIR)/lib -lBAT
GLIBS       += $(ROOTGLIBS) -L.-L$(BATINSTALLDIR)/lib -lBAT

CXSRCS      = TestSuite.cxx \
              ReleaseTestSuite.cxx \
              PerfTest.cxx \
              PerfTestVarPar.cxx \
              PerfSubTest.cxx \
              PerfTestMCMC.cxx \
              PerfTest1DFunction.cxx \
              PerfTest2DFunction.cxx

LIBA         = libBATTS.a
LIBSO        = libBATTS.so

CXXSRCS      = $(patsubst %.cxx,src/%.cxx,$(CXSRCS))

DICTSRCS     = $(patsubst %.cxx,src/%Dict.cxx,$(CXSRCS))

DICTHDRS     = $(patsubst %.cxx,src/%Dict.h,$(CXSRCS))

CXXOBJS      = $(patsubst %.cxx,obj/%.o,$(CXSRCS))

DICTOBJS     = $(patsubst %.cxx,obj/%Dict.o,$(CXSRCS))

EXEOBJS      =

GARBAGE      = $(CXXOBJS) $(EXEOBJS)  src/*Dict.*  obj/*.o  $(LIBA) $(LIBSO) link.d runTestSuite.exe 

all : $(LIBSO) runTestSuite

libs : $(LIBA) $(LIBSO)

link.d : $(patsubst %.cxx,include/%.h,$(CXSRCS))
	$(CXX) -MM $(CXXFLAGS) $(CXXSRCS) > link.d;

-include link.d

obj/%.o : src/%.cxx
	$(CXX) $(CXXFLAGS) -c $< -o $@

runTestSuite : runTestSuite.cxx 
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(LIBS) -lBATTS -c runTestSuite.cxx -o obj/runTestSuite.o
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(LIBS) -lBATTS  obj/runTestSuite.o -o runTestSuite.exe

$(LIBSO) : $(CXXOBJS) 
	$(CXX) $(SOFLAGS) $(LDFLAGS) $^ -o $(LIBSO)

$(LIBA) : $(CXXOBJS)
	$(AR)  $(ARFLAGS) $(LIBA) $(CXXOBJS)
	$(RANLIB) $(LIBA)

clean :
	$(RM) $(GARBAGE)

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

web :
	tar -cf web.tar \
	*.html \
	*.png \
	*.eps \
	*.ps \
	*.log
	gzip web.tar