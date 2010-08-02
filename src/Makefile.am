lib_LTLIBRARIES= libBAT.la
libBAT_la_LDFLAGS= -version-info 2:0:0

AM_CXXFLAGS= -Wall -fPIC  -I.. $(ROOTCFLAGS)
AM_LDFLAGS= -O
AM_SOFLAGS= -shared

if CUBA
AM_CXXFLAGS += $(CUBAINCLUDE)
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
library_include_HEADERS=../BAT/*.h
libBAT_la_SOURCES= \
   BCDataPoint.cxx \
   BCH2D.cxx \
   BCModel.cxx \
   BCParameter.cxx \
   BCDataSet.cxx \
   BCIntegrate.cxx \
   BCModelManager.cxx \
   BCEngineMCMC.cxx \
   BCLog.cxx \
   BCModelOutput.cxx \
   BCH1D.cxx \
   BCMath.cxx \
   BCGoFTest.cxx \
   BCAux.cxx \
   BCSummaryTool.cxx \
   BCSummaryPriorModel.cxx \
   \
   BCDataPointDict.cxx \
   BCH2DDict.cxx \
   BCModelDict.cxx \
   BCParameterDict.cxx \
   BCDataSetDict.cxx \
   BCIntegrateDict.cxx \
   BCModelManagerDict.cxx \
   BCEngineMCMCDict.cxx \
   BCLogDict.cxx \
   BCModelOutputDict.cxx \
   BCH1DDict.cxx \
   BCMathDict.cxx \
   BCGoFTestDict.cxx \
   BCAuxDict.cxx \
   BCSummaryToolDict.cxx \
   BCSummaryPriorModelDict.cxx

clean-local:
	rm -f *.o *.lo *Dict.*

%.o : %.cxx
	$(CXX) $(AM_CXXFLAGS) $(CXXFLAGS) -c $<

%Dict.cxx : ../BAT/%.h
	$(ROOTCINT) -f $@ -c  -I.. $<

install-exec-hook:
	cp libBAT.rootmap $(DESTDIR)$(libdir)
	$(CREATE_SO_LINK)
