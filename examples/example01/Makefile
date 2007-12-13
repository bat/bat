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

CXSRCS      = BCModelPol1.cxx 

CXXSRCS      = $(patsubst %.cxx,src/%.cxx,$(CXSRCS))

CXXOBJS      = $(patsubst %.cxx,%.o,$(CXSRCS))

EXEOBJS       = 

GARBAGE      = $(CXXOBJS) $(EXEOBJS) *.o *~ link.d runBATexample01.*.exe 

all : runBATexample01.1  runBATexample01.2  runBATexample01.3  runBATexample01.4  runBATexample01.5  runBATexample01.6  runBATexample01.7 

link.d : $(patsubst %.cxx,src/%.h,$(CXSRCS))
	$(CXX) -MM $(CXXFLAGS) $(CXXSRCS) > link.d; 

include link.d

%.o : src/%.cxx 
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean :
	$(RM) $(GARBAGE)

runBATexample01.1 : runBATexample01.1.c  $(CXXOBJS)
	$(CXX) $(CXXFLAGS) -c $<
	$(CXX) $(LDFLAGS) $(LIBS) runBATexample01.1.o $(CXXOBJS) -o runBATexample01.1.exe

runBATexample01.2 : runBATexample01.2.c  $(CXXOBJS)
	$(CXX) $(CXXFLAGS) -c $<
	$(CXX) $(LDFLAGS) $(LIBS) runBATexample01.2.o $(CXXOBJS) -o runBATexample01.2.exe

runBATexample01.3 : runBATexample01.3.c  $(CXXOBJS)
	$(CXX) $(CXXFLAGS) -c $<
	$(CXX) $(LDFLAGS) $(LIBS) runBATexample01.3.o $(CXXOBJS) -o runBATexample01.3.exe

runBATexample01.4 : runBATexample01.4.c  $(CXXOBJS)
	$(CXX) $(CXXFLAGS) -c $<
	$(CXX) $(LDFLAGS) $(LIBS) runBATexample01.4.o $(CXXOBJS) -o runBATexample01.4.exe

runBATexample01.5 : runBATexample01.5.c  $(CXXOBJS)
	$(CXX) $(CXXFLAGS) -c $<
	$(CXX) $(LDFLAGS) $(LIBS) runBATexample01.5.o $(CXXOBJS) -o runBATexample01.5.exe

runBATexample01.6 : runBATexample01.6.c  $(CXXOBJS)
	$(CXX) $(CXXFLAGS) -c $<
	$(CXX) $(LDFLAGS) $(LIBS) runBATexample01.6.o $(CXXOBJS) -o runBATexample01.6.exe

runBATexample01.7 : runBATexample01.7.c  $(CXXOBJS)
	$(CXX) $(CXXFLAGS) -c $<
	$(CXX) $(LDFLAGS) $(LIBS) runBATexample01.7.o $(CXXOBJS) -o runBATexample01.7.exe

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


