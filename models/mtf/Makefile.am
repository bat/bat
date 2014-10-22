CLEANFILES = *_rdict.* *_rdict_* *.rootmap

AM_CXXFLAGS= $(ROOTAUXCFLAGS) -Wall -fPIC -I$(top_srcdir) $(ROOTCFLAGS)

AM_LDFLAGS= -O
AM_SOFLAGS= -shared

if DEBUG
AM_CXXFLAGS += -g -O0 -Wall -Wextra -pedantic -Wno-long-long -Wno-unknown-pragmas
else
AM_CXXFLAGS += -O2 -Wno-unknown-pragmas
endif

if IS_DARWIN
CREATE_SO_LINK = cd $(DESTDIR)$(libdir); rm -f libBATmtf.so; ln -s libBATmtf.dylib libBATmtf.so
endif

lib_LTLIBRARIES= libBATmtf.la
libBATmtf_la_LDFLAGS= -version-info 2:0:0 -L$(ROOTLIBDIR)
libBATmtf_la_LIBADD = $(top_builddir)/src/libBAT.la $(ROOTLIBS)

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

nodist_libBATmtf_la_SOURCES= \
	libBATmtf_rdict.cxx

BUILT_SOURCES = libBATmtf_rdict.cxx

EXTRA_DIST = LinkDef.h

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
rootmap_DATA = .libs/libBATmtf.rootmap


if WITH_CLING
### ROOT-6: ###

rootmap_DATA += .libs/libBATmtf_rdict.pcm

libBATmtf_rdict.cxx: $(library_include_HEADERS) LinkDef.h
	$(ROOTCLING) -f $@.tmp -s libBATmtf@SHLIBEXT@ -rml libBATmtf@SHLIBEXT@ -rmf libBATmtf.rootmap -c $(CPPFLAGS) $(CXXFLAGS) -I$(includedir) $+
	@# Some magic to prefix header names with "$(PACKAGE)/", and only that, in dictionary and rootmap:
	$(GREP) -F -v '"'"`pwd`"'/",' $@.tmp | $(SED) 's|"\([^"]*/\)\?\([^/"]*[.]h\)",|"'$(PACKAGE)/'\2",| ; s|\\"\([^"]*/\)\?\([^/"]*[.]h\)\\"\\n"|\\"'$(PACKAGE)/'\2\\"\\n"|' > $@ && $(RM) $@.tmp
	$(SED) 's|\$$clingAutoload\$$\([^/""]\+/\)*|$$clingAutoload$$'$(PACKAGE)'/|; /.*DICTPAYLOAD(.*/,/.*)DICTPAYLOAD.*/ s|#include "\([^/"]\+/\)*\(.*\)"|#include <'$(PACKAGE)'/\2>|' -i $@
	$(SED) 's|\(header \+\)\([^ ].*/\)\?\([^ ].*[.]h\)|\1'$(PACKAGE)/'\3|' -i libBATmtf.rootmap

.libs/libBATmtf.rootmap: libBATmtf.la
	cp libBATmtf.rootmap libBATmtf_rdict.pcm .libs/

.libs/libBATmtf_rdict.pcm: .libs/libBATmtf.rootmap

else # not WITH_CLING
### ROOT-5: ###

libBATmtf_rdict.cxx: $(library_include_HEADERS) LinkDef.h
	$(ROOTCINT) -f $@ -c -p $(CPPFLAGS) $(CXXFLAGS) $+

.libs/libBATmtf.rootmap: libBATmtf.la | $(library_include_HEADERS) LinkDef.h
	$(RLIBMAP) -o $@ -l .libs/libBATmtf@SHLIBEXT@ -c $|

endif # WITH_CLING
