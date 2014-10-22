CLEANFILES = *_rdict.* *_rdict_* *.rootmap

lib_LTLIBRARIES = libBAT.la
libBAT_la_LDFLAGS = -version-info 6:0:0

# this fixes the linker error on ubuntu with gcc > 4.4 with missing reference to TCanvas(bool) etc.
libBAT_la_LDFLAGS += -L$(ROOTLIBDIR)
libBAT_la_LIBADD = $(ROOTLIBS)

AM_CXXFLAGS = $(ROOTAUXCFLAGS) -Wall -fPIC -I$(top_srcdir)/BAT $(ROOTCFLAGS)
AM_LDFLAGS  = -O
AM_SOFLAGS  = -shared

if CUBA
libBAT_la_LIBADD += -lcuba
endif

if DEBUG
AM_CXXFLAGS += -g -O0 -Wall -Wextra -pedantic -Wno-long-long -Wno-unknown-pragmas
else
AM_CXXFLAGS += -O2 -Wno-unknown-pragmas
endif

if IS_DARWIN
CREATE_SO_LINK = cd $(DESTDIR)$(libdir); rm -f libBAT.so; ln -s libBAT.dylib libBAT.so
endif

EXTRA_DIST = LinkDef.h

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
	$(patsubst %.h, $(top_srcdir)/BAT/%.h, $(incs))
libBAT_la_SOURCES= \
	$(incs:.h=.cxx)

nodist_libBAT_la_SOURCES= \
	libBAT_rdict.cxx

BUILT_SOURCES = libBAT_rdict.cxx

install-exec-hook:
	$(CREATE_SO_LINK)

# turn #include "file.h" into #include <file.h>
# work in tmp directory to allow VPATH build with read-only src directory
install-data-hook:
	mkdir -p src-tmp; \
	for i in $(DESTDIR)$(includedir)/BAT/*.h; do \
		sed -e 's/include \"\(BC.*.h\)\"/include <BAT\/\1>/' $$i > src-tmp/`basename $$i`.h; \
		if ! diff -q src-tmp/`basename $$i`.h "$$i" > /dev/null; then \
			echo "updating $$i"; \
			cp -f src-tmp/`basename $$i`.h "$$i"; \
			chmod 664 "$$i"; \
		fi; \
	done; \
	rm -rf src-tmp;


rootmapdir = $(libdir)
rootmap_DATA = .libs/libBAT.rootmap


if WITH_CLING
### ROOT-6: ###

rootmap_DATA += .libs/libBAT_rdict.pcm

libBAT_rdict.cxx: $(library_include_HEADERS) LinkDef.h
	$(ROOTCLING) -f $@.tmp -s libBAT@SHLIBEXT@ -rml libBAT@SHLIBEXT@ -rmf libBAT.rootmap -c $(CPPFLAGS) $(CXXFLAGS) -I$(includedir) $+
	@# Some magic to prefix header names with "$(PACKAGE)/", and only that, in dictionary and rootmap:
	$(GREP) -F -v '"'"`pwd`"'/",' $@.tmp | $(SED) 's|"\([^"]*/\)\?\([^/"]*[.]h\)",|"'$(PACKAGE)/'\2",| ; s|\\"\([^"]*/\)\?\([^/"]*[.]h\)\\"\\n"|\\"'$(PACKAGE)/'\2\\"\\n"|' > $@ && $(RM) $@.tmp
	$(SED) 's|\$$clingAutoload\$$\([^/""]\+/\)*|$$clingAutoload$$'$(PACKAGE)'/|; /.*DICTPAYLOAD(.*/,/.*)DICTPAYLOAD.*/ s|#include "\([^/"]\+/\)*\(.*\)"|#include <'$(PACKAGE)'/\2>|' -i $@
	$(SED) 's|\(header \+\)\([^ ].*/\)\?\([^ ].*[.]h\)|\1'$(PACKAGE)/'\3|' -i libBAT.rootmap

.libs/libBAT.rootmap: libBAT.la
	cp libBAT.rootmap libBAT_rdict.pcm .libs/

.libs/libBAT_rdict.pcm: .libs/libBAT.rootmap

else # not WITH_CLING
### ROOT-5: ###

libBAT_rdict.cxx: $(library_include_HEADERS) LinkDef.h
	$(ROOTCINT) -f $@ -c -p $(CPPFLAGS) $(CXXFLAGS) $+

.libs/libBAT.rootmap: libBAT.la | $(library_include_HEADERS) LinkDef.h
	$(RLIBMAP) -o $@ -l .libs/libBAT@SHLIBEXT@ -c $|

endif # WITH_CLING
