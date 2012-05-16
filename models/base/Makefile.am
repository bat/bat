lib_LTLIBRARIES= libBATmodels.la
libBATmodels_la_LDFLAGS= -version-info 3:0:0

AM_CXXFLAGS= -Wall -fPIC  -I../.. -I. -I../../BAT $(ROOTCFLAGS)

AM_LDFLAGS= -O
AM_SOFLAGS= -shared

if CUBA
#AM_CXXFLAGS += $(CUBAINCLUDE)
#AM_LDFLAGS += $(CUBALIB) -lcuba
#libBATmodels_la_LDFLAGS += $(CUBALIB) -lcuba
AM_LDFLAGS += -lcuba
libBATmodels_la_LDFLAGS += -lcuba
endif

if DEBUG
AM_CXXFLAGS += -g -O0 -Wextra -pedantic -Wno-long-long
else !DEBUG
AM_CXXFLAGS += -O2
endif !DEBUG

if IS_DARWIN
CREATE_SO_LINK = cd $(DESTDIR)$(libdir); rm -f libBATmodels.so; ln -s libBATmodels.dylib libBATmodels.so
endif

library_includedir= $(includedir)/BAT

incs = \
	BCEfficiencyFitter.h \
	BCGraphFitter.h \
	BCHistogramFitter.h \
	BCTemplateFitter.h \
	BCTemplateEnsembleTest.h

## RooStats specific options
if ROOSTATS
AM_LDFLAGS += -lRooFitCore -lRooFit -lRooStats -lFoam -lMathMore
incs += \
	BCRooInterface.h \
	BATCalculator.h
ADD_ROOINTERFACE_ROOTMAP = \
	rm -f libBATmodels.rootmap; \
	sed -e 's/$$/ libRooFit.so libRooFitCore.so libRooStats.so/' libBATmodels.rootmap.in > libBATmodels.rootmap; \
	cat BCRooInterface.rootmap >> libBATmodels.rootmap
endif
## End of RooStats specific options

library_include_HEADERS= $(incs)
libBATmodels_la_SOURCES= \
	$(incs:.h=.cxx)

#	$(incs:.h=Dict.cxx)

nodist_libBATmodels_la_SOURCES= \
	$(incs:.h=Dict.cxx)

BUILT_SOURCES = $(incs:.h=Dict.cxx)

clean-local:
	rm -f *.o *.lo *Dict.* libBATmodels.rootmap

%.o : %.cxx
	$(CXX) $(AM_CXXFLAGS) $(CXXFLAGS) -c $<

## RooStats specific targers
#if ROOSTATS
BATCalculatorDict.cxx : BATCalculator.h
	$(ROOTCINT) -f $@ -c  -I../.. -I. -I../../BAT $< BATCalculatorLinkDef.h
#endif
## End of RooStats specific targets

%Dict.cxx : %.h
	$(ROOTCINT) -f $@ -c  -I../.. -I. -I../../BAT $<

install-exec-hook:
	rm -rf libBATmodels.rootmap
	cp libBATmodels.rootmap.in libBATmodels.rootmap
	$(ADD_ROOINTERFACE_ROOTMAP)
	cp libBATmodels.rootmap $(DESTDIR)$(libdir)
	$(CREATE_SO_LINK)

install-data-hook:
	for i in $(incs); do \
		sed -e 's/include \"\(BC.*.h\)\"/include <BAT\/\1>/ ; s/include \"\.\.\/\.\.\/\(BAT\/BC.*.h\)\"/include <\1>/' $$i > $$i.h; \
		if ! diff -q "$$i.h" "$(DESTDIR)$(includedir)/BAT/$$i" > /dev/null; then \
			echo "updating $$i"; \
			cp -f "$$i.h" "$(DESTDIR)$(includedir)/BAT/$$i"; \
			chmod 664 "$(DESTDIR)$(includedir)/BAT/$$i"; \
		fi; \
		rm -f $$i.h; \
	done

uninstall-hook:
	rm -f $(DESTDIR)$(libdir)/libBATmodels.rootmap

