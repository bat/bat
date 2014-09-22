lib_LTLIBRARIES= libBATmtf.la
libBATmtf_la_LDFLAGS= -version-info 2:0:0

AM_CXXFLAGS= -Wall -fPIC -I$(top_srcdir) $(ROOTCFLAGS)

AM_LDFLAGS= -O
AM_SOFLAGS= -shared

if CUBA
AM_LDFLAGS += -lcuba
libBATmtf_la_LDFLAGS += -lcuba
endif

if DEBUG
AM_CXXFLAGS += -g -O0 -Wall -Wextra -pedantic -Wno-long-long -Wno-unknown-pragmas
else
AM_CXXFLAGS += -O2 -Wno-unknown-pragmas
endif

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
	$(ROOTCINT) -f $@ -c -I$(top_srcdir)/BAT $<

install-exec-hook:
	$(CREATE_SO_LINK)


install-data-hook:
	mkdir -p mtf-tmp; \
	for i in $(DESTDIR)$(includedir)/BAT/*.h; do \
		sed -e 's|include \"\(BC.*.h\)\"|include <BAT\/\1>| ; s/include \"\.\.\/\.\.\/\(BAT\/BC.*.h\)\"/include <\1>/' $$i > mtf-tmp/`basename $$i`.h; \
		if ! diff -q mtf-tmp/`basename $$i`.h "$$i" > /dev/null; then \
			echo "updating $$i"; \
			cp -f mtf-tmp/`basename $$i`.h "$$i"; \
			chmod 664 "$$i"; \
		fi; \
	done; \
	rm -rf mtf-tmp;

rootmapdir = $(libdir)
rootmap_DATA = libBATmtf.rootmap
