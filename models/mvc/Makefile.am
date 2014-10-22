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
CREATE_SO_LINK = cd $(DESTDIR)$(libdir); rm -f libBATmvc.so; ln -s libBATmvc.dylib libBATmvc.so
endif

lib_LTLIBRARIES= libBATmvc.la
libBATmvc_la_LDFLAGS= -version-info 0:0:0 -L$(ROOTLIBDIR)
libBATmvc_la_LIBADD = $(top_builddir)/src/libBAT.la $(ROOTLIBS)

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
	libBATmvc_rdict.cxx

BUILT_SOURCES = libBATmvc_rdict.cxx

EXTRA_DIST = LinkDef.h

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
rootmap_DATA = .libs/libBATmvc.rootmap


if WITH_CLING
### ROOT-6: ###

rootmap_DATA += .libs/libBATmvc_rdict.pcm

libBATmvc_rdict.cxx: $(library_include_HEADERS) LinkDef.h
	$(ROOTCLING) -f $@.tmp -s libBATmvc@SHLIBEXT@ -rml libBATmvc@SHLIBEXT@ -rmf libBATmvc.rootmap -c $(CPPFLAGS) $(CXXFLAGS) -I$(includedir) $+
	@# Some magic to prefix header names with "$(PACKAGE)/", and only that, in dictionary and rootmap:
	$(GREP) -F -v '"'"`pwd`"'/",' $@.tmp | $(SED) 's|"\([^"]*/\)\?\([^/"]*[.]h\)",|"'$(PACKAGE)/'\2",| ; s|\\"\([^"]*/\)\?\([^/"]*[.]h\)\\"\\n"|\\"'$(PACKAGE)/'\2\\"\\n"|' > $@ && $(RM) $@.tmp
	$(SED) 's|\$$clingAutoload\$$\([^/""]\+/\)*|$$clingAutoload$$'$(PACKAGE)'/|; /.*DICTPAYLOAD(.*/,/.*)DICTPAYLOAD.*/ s|#include "\([^/"]\+/\)*\(.*\)"|#include <'$(PACKAGE)'/\2>|' -i $@
	$(SED) 's|\(header \+\)\([^ ].*/\)\?\([^ ].*[.]h\)|\1'$(PACKAGE)/'\3|' -i libBATmvc.rootmap

.libs/libBATmvc.rootmap: libBATmvc.la
	cp libBATmvc.rootmap libBATmvc_rdict.pcm .libs/

.libs/libBATmvc_rdict.pcm: .libs/libBATmvc.rootmap

else # not WITH_CLING
### ROOT-5: ###

libBATmvc_rdict.cxx: $(library_include_HEADERS) LinkDef.h
	$(ROOTCINT) -f $@ -c -p $(CPPFLAGS) $(CXXFLAGS) $+

.libs/libBATmvc.rootmap: libBATmvc.la | $(library_include_HEADERS) LinkDef.h
	$(RLIBMAP) -o $@ -l .libs/libBATmvc@SHLIBEXT@ -c $|

endif # WITH_CLING
