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
CXXFLAGS    += -I./include -I./src -I./
LIBS        += $(ROOTLIBS) -L. -lbc -lcuba 
GLIBS       += $(ROOTGLIBS) -L. -lbc

CXSRCS      = BCParameter.cxx \
		BCModel.cxx \
		BCDataPoint.cxx \
		BCDataSet.cxx \
		BCH1D.cxx \
		BCH2D.cxx \
		BCIntegrate.cxx \
		BCModelManager.cxx \
		BCLog.cxx \

CXXSRCS      = $(patsubst %.cxx,src/%.cxx,$(CXSRCS))

CXXOBJS      = $(patsubst %.cxx,obj/%.o,$(CXSRCS))

EXEOBJS       = 

GARBAGE      = $(CXXOBJS) $(EXEOBJS)  libbc.so BayesianCalculator.o

all : libbc.so

link.d : $(patsubst %.cxx,include/%.h,$(CXSRCS))
	$(CXX) -MM $(CXXFLAGS) $(CXXSRCS) > link.d; 

include link.d

obj/%.o : src/%.cxx 
	$(CXX) $(CXXFLAGS) -c $< -o $@

libbc.so : $(CXXOBJS)  
	$(CXX) $(SOFLAGS) $(LDFLAGS) $^ -o $@

clean :
	$(RM) $(GARBAGE)

BaysianCalculator : BayesianCalculator.c libbc.so
	$(CXX) $(CXXFLAGS) -c $<
	$(CXX) $(LDFLAGS) $(ROOTLIBS) BayesianCalculator.o ${PWD}/libbc.so -o BayesianCalculator.exe 
	cp BayesianCalculator.exe exe/.

test : BayesianCalculator.c 
	$(CXX) $(CXXFLAGS) -c BayesianCalculator.c 
	$(CXX) $(LDFLAGS) $(ROOTLIBS) BayesianCalculator.o ${PWD}/libbc.so -o BayesianCalculator_2.exe 

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


