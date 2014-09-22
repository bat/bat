lib_LTLIBRARIES= libBATmvc.la
libBATmvc_la_LDFLAGS= -version-info 0:0:0

AM_CXXFLAGS= -Wall -fPIC -I$(top_srcdir) $(ROOTCFLAGS)

AM_LDFLAGS= -O
AM_SOFLAGS= -shared

if CUBA
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

library_include_HEADERS= $(incs)
libBATmvc_la_SOURCES= \
	$(incs:.h=.cxx)

nodist_libBATmvc_la_SOURCES= \
	$(incs:.h=Dict.cxx)

BUILT_SOURCES = $(incs:.h=Dict.cxx)

clean-local:
	rm -f *.o *.lo *Dict.*

%.o : %.cxx
	$(CXX) $(AM_CXXFLAGS) $(CXXFLAGS) -c $<

%Dict.cxx : %.h
	$(ROOTCINT) -f $@ -c -I$(top_srcdir) $<

install-exec-hook:
	$(CREATE_SO_LINK)

install-data-hook:
	mkdir -p mvc-tmp; \
	for i in $(DESTDIR)$(includedir)/BAT/*.h; do \
		sed -e 's|include \"\(BC.*.h\)\"|include <BAT\/\1>| ; s/include \"\.\.\/\.\.\/\(BAT\/BC.*.h\)\"/include <\1>/' $$i > mvc-tmp/`basename $$i`.h; \
		if ! diff -q mvc-tmp/`basename $$i`.h "$$i" > /dev/null; then \
			echo "updating $$i"; \
			cp -f mvc-tmp/`basename $$i`.h "$$i"; \
			chmod 664 "$$i"; \
		fi; \
	done; \
	rm -rf mvc-tmp;

rootmapdir = $(libdir)
rootmap_DATA = libBATmvc.rootmap
