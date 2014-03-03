lib_LTLIBRARIES= libBATmvc.la
libBATmvc_la_LDFLAGS= -version-info 0:0:0

AM_CXXFLAGS= -Wall -fPIC  -I../.. -I. -I../../BAT $(ROOTCFLAGS)

AM_LDFLAGS= -O
AM_SOFLAGS= -shared

if CUBA
#AM_CXXFLAGS += $(CUBAINCLUDE)
#AM_LDFLAGS += $(CUBALIB) -lcuba
#libBATmvc_la_LDFLAGS += $(CUBALIB) -lcuba
AM_LDFLAGS += -lcuba
libBATmvc_la_LDFLAGS += -lcuba
endif

if DEBUG
AM_CXXFLAGS += -g -O0 -Wall -Wextra -pedantic -Wno-long-long -Wno-unknown-pragmas
else !DEBUG
AM_CXXFLAGS += -O2 -Wno-unknown-pragmas
endif !DEBUG

if IS_DARWIN
CREATE_SO_LINK = cd $(DESTDIR)$(libdir); rm -f libBATmvc.so; ln -s libBATmvc.dylib libBATmvc.so
endif

library_includedir= $(includedir)/BAT

incs = \
	BCMVCombination.h \
	BCMVCUncertainty.h \
	BCMVCMeasurement.h \
	BCMVCObservable.h \
	BCMVCDataModel.h \
	BCMVCPhysicsModel.h

## RooStats specific options
if ROOSTATS
AM_LDFLAGS += -lRooFitCore -lRooFit -lRooStats -lFoam -lMathMore
endif
## End of RooStats specific options

library_include_HEADERS= $(incs)
libBATmvc_la_SOURCES= \
	$(incs:.h=.cxx)

#	$(incs:.h=Dict.cxx)

nodist_libBATmvc_la_SOURCES= \
	$(incs:.h=Dict.cxx)

BUILT_SOURCES = $(incs:.h=Dict.cxx)

clean-local:
	rm -f *.o *.lo *Dict.*

%.o : %.cxx
	$(CXX) $(AM_CXXFLAGS) $(CXXFLAGS) -c $<

%Dict.cxx : %.h
	$(ROOTCINT) -f $@ -c  -I../.. -I. -I../../BAT $<

install-exec-hook:
	cp libBATmvc.rootmap $(DESTDIR)$(libdir)
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
	rm -f $(DESTDIR)$(libdir)/libBATmvc.rootmap

