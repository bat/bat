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
CXXFLAGS    += -I.

########################################################
# If your version of BAT was not compiled with cuba and
# you don't have Cuba installed in your system you might
# need to remove '-lcuba' from the two following lines.
# Otherwise the example won't compile.
###
LIBS        += $(ROOTLIBS)  -lBAT -lBATmodels
# -lcuba
GLIBS       += $(ROOTGLIBS) -lBAT -lBATmodels
# -lcuba

CXXSRCS      = \
	Pol1Asymm.cxx \
	Pol2Asymm.cxx

CXXOBJS      = $(patsubst %.cxx,%.o,$(CXXSRCS))
EXEOBJS      =
GARBAGE      = $(CXXOBJS) $(EXEOBJS) *.o *~ link.d asymmErrorsFitterExample

all : asymmErrorsFitterExample

link.d : $(patsubst %.cxx,%.h,$(CXXSRCS))
	$(CXX) -MM $(CXXFLAGS) $(CXXSRCS) > link.d;

include link.d

%.o : %.cxx
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean :
	$(RM) $(GARBAGE)

asymmErrorsFitterExample : asymmErrorsFitterExample.cxx $(CXXOBJS)
	$(CXX) $(CXXFLAGS) -c $<
	$(CXX) $(LDFLAGS) $(LIBS) $@.o $(CXXOBJS) -o $@

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


