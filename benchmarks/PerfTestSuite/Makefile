# Programs
CXX          = g++
CXXFLAGS     = -g -O0 -Wall -Wextra -pedantic -Werror -fPIC
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
CXXFLAGS    += $(shell root-config --cflags)
CXXFLAGS    += -I./ -I./include -I/src $(shell bat-config --cflags)
LIBS        += -L. $(shell bat-config --libs)

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

EXEOBJS      = runTestSuite

GARBAGE      = $(CXXOBJS) $(EXEOBJS)  src/*Dict.*  obj/*.o  $(LIBA) $(LIBSO) link.d runTestSuite

RESULTS      = *.png *.pdf *.log *.tst *.html *.php *.htm

all : directory $(LIBSO) runTestSuite

libs : $(LIBA) $(LIBSO)

link.d : $(patsubst %.cxx,include/%.h,$(CXSRCS))
	$(CXX) -MM $(CXXFLAGS) $(CXXSRCS) > link.d;

-include link.d

obj/%.o : src/%.cxx
	$(CXX) $(CXXFLAGS) -c $< -o $@

runTestSuite : runTestSuite.cxx
	$(CXX) $(CXXFLAGS) $(LDFLAGS)  $(LIBS) -lBATTS -c runTestSuite.cxx -o obj/runTestSuite.o
	$(CXX) $(CXXFLAGS) $(LDFLAGS)  obj/runTestSuite.o -lBATTS $(LIBS)   -o runTestSuite

$(LIBSO) : $(CXXOBJS)
	$(CXX) $(SOFLAGS) $(LDFLAGS) $(LIBS) $^ -o $(LIBSO)

$(LIBA) : $(CXXOBJS)
	$(AR)  $(ARFLAGS) $(LIBS) $(LIBA) $(CXXOBJS)
	$(RANLIB) $(LIBA)

clean :
	$(RM) $(GARBAGE)

resclean :
	$(RM) $(RESULTS)

directory :
	mkdir -p obj

print :
	echo compiler  : $(CXX)
	echo compiler  : $(CXSRCS)
	echo c++ srcs  : $(CXXSRCS)
	echo c++ objs  : $(CXXOBJS)
	echo c++ flags : $(CXXFLAGS)
	echo libs      : $(LIBS)
	echo so flags  : $(SOFLAGS)

web :
	tar -zcf web.tar.gz \
	*.php \
	*.png \
	*.pdf \
	*.log
