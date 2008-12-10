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
LIBS        += $(ROOTLIBS) -L. -L../../ -lBAT -lcuba
GLIBS       += $(ROOTGLIBS) -L. -L../../ -lBAT -lcuba

CXSRCS       = 

CXXSRCS      = $(patsubst %.cxx,src/%.cxx,$(CXSRCS))

CXXOBJS      = $(patsubst %.cxx,%.o,$(CXSRCS))

EXEOBJS      = 

GARBAGE      = $(CXXOBJS) $(EXEOBJS) *.o *~ link.d runBATexample07.exe 

all : runBATexample07

link.d : $(patsubst %.cxx,src/%.h,$(CXSRCS))
	$(CXX) -MM $(CXXFLAGS) $(CXXSRCS) > link.d; 

#include link.d

%.o : src/%.cxx 
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean :
	$(RM) $(GARBAGE)

runBATexample07 : runBATexample07.c $(CXXOBJS)
	$(CXX) $(CXXFLAGS) -c $<
	$(CXX) $(LDFLAGS) $(LIBS) runBATexample07.o $(CXXOBJS) -o runBATexample07.exe

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


