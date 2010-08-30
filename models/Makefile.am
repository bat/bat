lib_LTLIBRARIES= libBATmodels.la
libBATmodels_la_LDFLAGS= -version-info 0:0:0

AM_CXXFLAGS= -Wall -fPIC  -I.. $(ROOTCFLAGS)

AM_LDFLAGS= -O
AM_SOFLAGS= -shared

if CUBA
AM_CXXFLAGS += $(CUBAINCLUDE)
AM_LDFLAGS += -lcuba
libBATmodels_la_LDFLAGS += -lcuba
endif

if DEBUG
AM_CXXFLAGS += -g -O0
else !DEBUG
AM_CXXFLAGS += -O2
endif !DEBUG

if IS_DARWIN
CREATE_SO_LINK = cd $(DESTDIR)$(libdir); rm -f libBATmodels.so; ln -s libBATmodels.dylib libBATmodels.so
endif

library_includedir= $(includedir)/BAT
library_include_HEADERS= \
	BCEfficiencyFitter.h \
	BCGraphFitter.h \
	BCHistogramFitter.h \
	BCTemplateFitter.h \
	BCTemplateEnsembleTest.h

libBATmodels_la_SOURCES= \
	BCEfficiencyFitter.cxx \
	BCGraphFitter.cxx \
	BCHistogramFitter.cxx \
	BCTemplateFitter.cxx \
	BCTemplateEnsembleTest.cxx \
	\
	BCEfficiencyFitterDict.cxx \
	BCGraphFitterDict.cxx \
	BCHistogramFitterDict.cxx \
	BCTemplateFitterDict.cxx \
	BCTemplateEnsembleTestDict.cxx

## RooStats specific options
if ROOSTATS
AM_LDFLAGS += -lRooFitCore -lRooFit -lRooStats -lFoam
libBATmodels_la_SOURCES += \
	BCRooInterface.cxx \
	BATCalculator.cxx \
	\
	BCRooInterfaceDict.cxx \
	BATCalculatorDict.cxx
library_include_HEADERS += \
	BCRooInterface.h \
	BATCalculator.h
ADD_ROOINTERFACE_ROOTMAP = \
	rm -f libBATmodels.rootmap; \
	sed -e 's/$$/ libRooFit.so libRooFitCore.so libRooStats.so/' libBATmodels.rootmap.in > libBATmodels.rootmap; \
	cat BCRooInterface.rootmap >> libBATmodels.rootmap
endif
## End of RooStats specific options

clean-local:
	rm -f *.o *.lo *Dict.* libBATmodels.rootmap

%.o : %.cxx
	$(CXX) $(AM_CXXFLAGS) $(CXXFLAGS) -c $<

## RooStats specific targers
#if ROOSTATS
BATCalculatorDict.cxx : BATCalculator.h
	$(ROOTCINT) -f $@ -c  -I.. $< BATCalculatorLinkDef.h
#endif
## End of RooStats specific targets

%Dict.cxx : %.h
	$(ROOTCINT) -f $@ -c  -I.. $<

install-exec-hook:
	rm -rf libBATmodels.rootmap
	cp libBATmodels.rootmap.in libBATmodels.rootmap
	$(ADD_ROOINTERFACE_ROOTMAP)
	cp libBATmodels.rootmap $(DESTDIR)$(libdir)
	
	cp $(library_include_HEADERS) $(DESTDIR)$(includedir)/BAT
	$(CREATE_SO_LINK)


