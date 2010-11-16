lib_LTLIBRARIES= libBAT.la
libBAT_la_LDFLAGS= -version-info 3:0:0

AM_CXXFLAGS= -Wall -fPIC  -I.. $(ROOTCFLAGS)
AM_LDFLAGS= -O
AM_SOFLAGS= -shared

if CUBA
#AM_CXXFLAGS += $(CUBAINCLUDE)
#AM_LDFLAGS += $(CUBALIB) -lcuba
#libBAT_la_LDFLAGS += $(CUBALIB) -lcuba
AM_LDFLAGS += -lcuba
libBAT_la_LDFLAGS += -lcuba
endif

if DEBUG
AM_CXXFLAGS += -g -O0
else !DEBUG
AM_CXXFLAGS += -O2
endif !DEBUG

if IS_DARWIN
CREATE_SO_LINK= cd $(DESTDIR)$(libdir); rm -f libBAT.so; ln -s libBAT.dylib libBAT.so
endif

library_includedir= $(includedir)/BAT

incs= \
	BCDataPoint.h \
	BCDataSet.h \
	BCH1D.h \
	BCH2D.h \
	BCParameter.h \
	BCEngineMCMC.h \
	BCIntegrate.h \
	BCModel.h \
	BCModelManager.h \
	BCModelOutput.h \
	BCGoFTest.h \
	BCLog.h \
	BCMath.h \
	BCAux.h \
	BCSummaryTool.h \
	BCSummaryPriorModel.h

library_include_HEADERS= \
	$(patsubst %.h, ../BAT/%.h, $(incs)) \
	../BAT/BCErrorCodes.h
libBAT_la_SOURCES= \
	$(incs:.h=.cxx)

nodist_libBAT_la_SOURCES= \
	$(incs:.h=Dict.cxx)

BUILT_SOURCES = $(incs:.h=Dict.cxx)

clean-local:
	rm -f *.o *.lo *Dict.*

%.o : %.cxx
	$(CXX) $(AM_CXXFLAGS) $(CXXFLAGS) -c $<

%Dict.cxx : ../BAT/%.h
	$(ROOTCINT) -f $@ -c  -I.. $<

install-exec-hook:
	cp libBAT.rootmap $(DESTDIR)$(libdir)
	$(CREATE_SO_LINK)
