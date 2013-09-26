lib_LTLIBRARIES= libBAT.la
libBAT_la_LDFLAGS= -version-info 6:0:0

# this fixes the linker error on ubuntu with gcc > 4.4 with missing reference to TCanvas(bool) etc.
libBAT_la_LDFLAGS += -L$(ROOTLIBDIR)
libBAT_la_LIBADD = $(ROOTLIBS)

AM_CXXFLAGS= -Wall -fPIC  -I.. -I. -I../BAT $(ROOTCFLAGS)
AM_LDFLAGS= -O
AM_SOFLAGS= -shared

if CUBA
libBAT_la_LIBADD += -lcuba
endif

if DEBUG
AM_CXXFLAGS += -g -O0 -Wall -Wextra -pedantic -Wno-long-long -Wno-unknown-pragmas
else
AM_CXXFLAGS += -O2 -Wno-unknown-pragmas
endif

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
	BCParameterSet.h \
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
	$(patsubst %.h, ../BAT/%.h, $(incs))
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

install-data-hook:
	for i in $(library_include_HEADERS); do \
		sed -e 's/include \"\(BC.*.h\)\"/include <BAT\/\1>/' $$i > $$i.h; \
		if ! diff -q "$$i.h" "$(DESTDIR)$(includedir)/BAT/$$i" > /dev/null; then \
			echo "updating $$i"; \
			cp -f "$$i.h" "$(DESTDIR)$(includedir)/BAT/$$i"; \
			chmod 664 "$(DESTDIR)$(includedir)/BAT/$$i"; \
		fi; \
		rm -f $$i.h; \
	done

uninstall-hook:
	rm -f $(DESTDIR)$(libdir)/libBAT.rootmap
