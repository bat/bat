lib_LTLIBRARIES= libBATmtf.la
libBATmtf_la_LDFLAGS= -version-info 1:0:1

AM_CXXFLAGS= -Wall -fPIC  -I../.. -I. -I../../BAT $(ROOTCFLAGS)

AM_LDFLAGS= -O
AM_SOFLAGS= -shared

if CUBA
#AM_CXXFLAGS += $(CUBAINCLUDE)
#AM_LDFLAGS += $(CUBALIB) -lcuba
#libBATmtf_la_LDFLAGS += $(CUBALIB) -lcuba
AM_LDFLAGS += -lcuba
libBATmtf_la_LDFLAGS += -lcuba
endif

if DEBUG
AM_CXXFLAGS += -g -O0 -Wextra -pedantic -Wno-long-long
else !DEBUG
AM_CXXFLAGS += -O2
endif !DEBUG

if IS_DARWIN
CREATE_SO_LINK = cd $(DESTDIR)$(libdir); rm -f libBATmtf.so; ln -s libBATmtf.dylib libBATmtf.so
endif

library_includedir= $(includedir)/BAT

incs = \
	BCMTFAnalysisFacility.h \
	BCMTFChannel.h \
	BCMTFComparisonTool.h \
	BCMTF.h \
	BCMTFProcess.h \
	BCMTFSystematic.h \
	BCMTFSystematicVariation.h \
	BCMTFTemplate.h

## RooStats specific options
if ROOSTATS
AM_LDFLAGS += -lRooFitCore -lRooFit -lRooStats -lFoam -lMathMore
endif
## End of RooStats specific options

library_include_HEADERS= $(incs)
libBATmtf_la_SOURCES= \
	$(incs:.h=.cxx)

#	$(incs:.h=Dict.cxx)

nodist_libBATmtf_la_SOURCES= \
	$(incs:.h=Dict.cxx)

BUILT_SOURCES = $(incs:.h=Dict.cxx)

clean-local:
	rm -f *.o *.lo *Dict.*

%.o : %.cxx
	$(CXX) $(AM_CXXFLAGS) $(CXXFLAGS) -c $<

%Dict.cxx : %.h
	$(ROOTCINT) -f $@ -c  -I../.. -I. -I../../BAT $<

install-exec-hook:
	cp libBATmtf.rootmap $(DESTDIR)$(libdir)
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
	rm -f $(DESTDIR)$(libdir)/libBATmtf.rootmap

