CLEANFILES = *_rdict.* *_rdict_* *.rootmap LinkDef.h

AM_CXXFLAGS= $(ROOTAUXCFLAGS) -Wall -fPIC -I$(top_srcdir) $(ROOTCFLAGS)
AM_LDFLAGS= -O
AM_SOFLAGS= -shared

lib_LTLIBRARIES= libBATmodels.la
libBATmodels_la_LDFLAGS= -version-info 4:0:0 -L$(ROOTLIBDIR)
libBATmodels_la_LIBADD = $(top_builddir)/src/libBAT.la $(ROOTLIBS)

if DEBUG
AM_CXXFLAGS += -g -O0 -Wall -Wextra -pedantic -Wno-long-long -Wno-unknown-pragmas
else
AM_CXXFLAGS += -O2 -Wno-unknown-pragmas
endif

if IS_DARWIN
CREATE_SO_LINK = cd $(DESTDIR)$(libdir); rm -f libBATmodels.so; ln -s libBATmodels.dylib libBATmodels.so
endif

library_includedir= $(includedir)/BAT

incs = \
	BCFitter.h \
	BCEfficiencyFitter.h \
	BCGraphFitter.h \
	BCHistogramFitter.h

linkdefs = LDef.h

## RooStats specific options
if ROOSTATS
incs += \
	BCRooInterface.h \
	BATCalculator.h
linkdefs += LDef_roostats.h
endif
## End of RooStats specific options

LinkDef.h: $(linkdefs)
	cat $+ > $@

library_include_HEADERS= $(incs)
libBATmodels_la_SOURCES= \
	$(incs:.h=.cxx)

nodist_libBATmodels_la_SOURCES= \
	libBATmodels_rdict.cxx

BUILT_SOURCES = libBATmodels_rdict.cxx LinkDef.h

EXTRA_DIST = LDef.h LDef_roostats.h

install-exec-hook:
	$(CREATE_SO_LINK)

install-data-hook:
	mkdir -p base-tmp; \
	for i in $(DESTDIR)$(includedir)/BAT/*.h; do \
		sed -e 's|include \"\(BC.*.h\)\"|include <BAT\/\1>| ; s/include \"\.\.\/\.\.\/\(BAT\/BC.*.h\)\"/include <\1>/' $$i > base-tmp/`basename $$i`.h; \
		if ! diff -q base-tmp/`basename $$i`.h "$$i" > /dev/null; then \
			echo "updating $$i"; \
			cp -f base-tmp/`basename $$i`.h "$$i"; \
			chmod 664 "$$i"; \
		fi; \
	done; \
	rm -rf base-tmp;

rootmapdir = $(libdir)
rootmap_DATA = .libs/libBATmodels.rootmap

if WITH_CLING
### ROOT-6: ###

rootmap_DATA += .libs/libBATmodels_rdict.pcm

libBATmodels_rdict.cxx: $(library_include_HEADERS) LinkDef.h
	$(ROOTCLING) -f $@.tmp -s libBATmodels@SHLIBEXT@ -rml libBATmodels@SHLIBEXT@ -rmf libBATmodels.rootmap -c $(CPPFLAGS) $(CXXFLAGS) -I$(includedir) $+
	@# Some magic to prefix header names with "$(PACKAGE)/", and only that, in dictionary and rootmap:
	$(GREP) -F -v '"'"`pwd`"'/",' $@.tmp | $(SED) 's|"\([^"]*/\)\?\([^/"]*[.]h\)",|"'$(PACKAGE)/'\2",| ; s|\\"\([^"]*/\)\?\([^/"]*[.]h\)\\"\\n"|\\"'$(PACKAGE)/'\2\\"\\n"|' > $@ && $(RM) $@.tmp
	$(SED) 's|\$$clingAutoload\$$\([^/""]\+/\)*|$$clingAutoload$$'$(PACKAGE)'/|; /.*DICTPAYLOAD(.*/,/.*)DICTPAYLOAD.*/ s|#include "\([^/"]\+/\)*\(.*\)"|#include <'$(PACKAGE)'/\2>|' -i $@
	$(SED) 's|\(header \+\)\([^ ].*/\)\?\([^ ].*[.]h\)|\1'$(PACKAGE)/'\3|' -i libBATmodels.rootmap

.libs/libBATmodels.rootmap: libBATmodels.la
	cp libBATmodels.rootmap libBATmodels_rdict.pcm .libs/

.libs/libBATmodels_rdict.pcm: .libs/libBATmodels.rootmap

else # not WITH_CLING
### ROOT-5: ###

libBATmodels_rdict.cxx: $(library_include_HEADERS) LinkDef.h
	$(ROOTCINT) -f $@ -c -p $(CPPFLAGS) $(CXXFLAGS) $+

.libs/libBATmodels.rootmap: libBATmodels.la | $(library_include_HEADERS) LinkDef.h
	$(RLIBMAP) -o $@ -l .libs/libBATmodels@SHLIBEXT@ -c $|

endif # WITH_CLING
