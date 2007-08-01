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
CXXFLAGS    += -I./src/ -I./ -I../../include/ 
LIBS        += $(ROOTLIBS) -L. -lbc -lcuba
GLIBS       += $(ROOTGLIBS) -L. -lbc -lcuba

CXSRCS      = BCModelPol0.cxx \
		BCModelPol1.cxx \
		BCModelPol2.cxx 

CXXSRCS      = $(patsubst %.cxx,src/%.cxx,$(CXSRCS))

CXXOBJS      = $(patsubst %.cxx,%.o,$(CXSRCS))

EXEOBJS       = 

GARBAGE      = $(CXXOBJS) $(EXEOBJS) libex2.so BayesianCalculator.o

all : BaysianCalculator

link.d : $(patsubst %.cxx,src/%.h,$(CXSRCS))
	$(CXX) -MM $(CXXFLAGS) $(CXXSRCS) > link.d; 

include link.d

%.o : src/%.cxx 
	$(CXX) $(CXXFLAGS) -c $< -o $@

libex2.so : $(CXXOBJS)  
	$(CXX) $(SOFLAGS) $(LDFLAGS) $^ -o $@

clean :
	$(RM) $(GARBAGE)

BaysianCalculator : BayesianCalculator.c libex2.so
	$(CXX) $(CXXFLAGS) -c $<
	$(CXX) $(LDFLAGS) $(ROOTLIBS) -L../../ -lbc BayesianCalculator.o -lcuba ${PWD}/libex2.so -o BayesianCalculator.exe 

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


