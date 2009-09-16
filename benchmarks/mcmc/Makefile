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
CXXFLAGS    += -I./src/ -I./ -I$(BATINSTALLDIR)/include
LIBS        += $(ROOTLIBS) -L. -L$(BATINSTALLDIR)/lib -lBAT
GLIBS       += $(ROOTGLIBS) -L. -L$(BATINSTALLDIR)/lib -lBAT

CXSRCS      = BCBenchmarkMCMC3D.cxx BCBenchmarkMCMC2D.cxx BCBenchmarkMCMC.cxx

CXXSRCS      = $(patsubst %.cxx,src/%.cxx,$(CXSRCS))

CXXOBJS      = $(patsubst %.cxx,%.o,$(CXSRCS))

EXEOBJS       =

GARBAGE      = $(CXXOBJS) $(EXEOBJS) *.o *.log link.d benchmark

all : benchmark

link.d : $(patsubst %.cxx,src/%.h,$(CXSRCS))
	$(CXX) -MM $(CXXFLAGS) $(CXXSRCS) > link.d;

include link.d

%.o : src/%.cxx
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean :
	$(RM) $(GARBAGE)

benchmark : benchmark.cxx $(CXXOBJS)
	$(CXX) $(CXXFLAGS) -c $<
	$(CXX) $(LDFLAGS) $(LIBS) $@.o $(CXXOBJS) -o $@

print :
	@echo compiler  : $(CXX)
	@echo compiler  : $(CXSRCS)
	@echo c++ srcs  : $(CXXSRCS)
	@echo c++ objs  : $(CXXOBJS)
	@echo c++ flags : $(CXXFLAGS)
	@echo libs      : $(LIBS)
	@echo so flags  : $(SOFLAGS)
	@echo rootlibs  : $(ROOTLIBS)
	@echo rootglibs : $(ROOTGLIBS)

