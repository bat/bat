CLEANFILES = *_rdict.* *_rdict_* *.rootmap

clean-local:
	-rm -rf .includes

lib_LTLIBRARIES= libBATmvc.la
libBATmvc_la_LDFLAGS= -version-info 0:0:0 -L$(ROOTLIBDIR)
libBATmvc_la_LIBADD = $(top_builddir)/src/libBAT.la $(ROOTLIBS)

AM_CXXFLAGS = $(ROOTAUXCFLAGS) -Wall -fPIC -I.includes -I$(top_builddir)/src/.includes $(ROOTCFLAGS)
AM_LDFLAGS  = -O
AM_SOFLAGS  = -shared

if DEBUG
AM_CXXFLAGS += -g -O0 -Wall -Wextra -pedantic -Wno-long-long -Wno-unknown-pragmas
else
AM_CXXFLAGS += -O2 -Wno-unknown-pragmas
endif

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

library_include_HEADERS= $(incs:%.h=.includes/$(PACKAGE)/%.h)

$(library_include_HEADERS): .includes/$(PACKAGE)/%.h: %.h
	mkdir -p .includes/$(PACKAGE)
	sed 's|include \"\(BC.*.h\)\"|include <'$(PACKAGE)'\/\1>| ; s/include \"\.\.\/\.\.\/\('$(PACKAGE)'\/BC.*.h\)\"/include <\1>/' $< > $@

libBATmvc_la_SOURCES= \
	$(incs:.h=.cxx)

nodist_libBATmvc_la_SOURCES= \
	libBATmvc_rdict.cxx

BUILT_SOURCES = libBATmvc_rdict.cxx

EXTRA_DIST = LinkDef.h

install-exec-hook:
	$(CREATE_SO_LINK)

rootmapdir = $(libdir)
rootmap_DATA = .libs/libBATmvc.rootmap


if WITH_CLING
### ROOT-6: ###

rootmap_DATA += .libs/libBATmvc_rdict.pcm

# short hand to make rule more concise, rtmap is not a standard automake variable
rtmap = libBATmvc.rootmap
libBATmvc_rdict.cxx: $(library_include_HEADERS) LinkDef.h
	$(ROOTCLING) -f $@.tmp -s libBATmvc@SHLIBEXT@ -rml libBATmvc@SHLIBEXT@ -rmf libBATmvc.rootmap -c $(CPPFLAGS) $(CXXFLAGS) -I.includes -I$(top_builddir)/src/.includes -I$(includedir) $+
	@# Some magic to prefix header names with "$(PACKAGE)/", and only that, in dictionary and rootmap:
	@# avoid sed -i: not portable between Mac OS X and linux => use back-up files
	$(GREP) -F -v '"'"`pwd`"'/",' $@.tmp | grep -v '^\"\([^/"]*/\)*\.includes\",*$$' | $(SED) 's|"\([^/"]*/\)*\([^/"]*[.]h\)",|"'$(PACKAGE)/'\2",| ; s|\\"\([^/"]*/\)*\([^/"]*[.]h\)\\"\\n"|\\"'$(PACKAGE)/'\2\\"\\n"|' > $@.tmp2
	$(SED) 's|\$$clingAutoload\$$\([^/"]*/\)*|$$clingAutoload$$'$(PACKAGE)'/|; /.*DICTPAYLOAD(.*/,/.*)DICTPAYLOAD.*/ s|#include "\([^/"]*/\)*\(.*\)"|#include <'$(PACKAGE)'/\2>|' $@.tmp2 > $@ && $(RM) $@.tmp $@.tmp2
	$(SED) 's|\(header  *\)\([^ ].*/\)*\([^ ].*[.]h\)|\1'$(PACKAGE)/'\3|' $(rtmap) > $(rtmap).tmp && mv $(rtmap).tmp $(rtmap)

.libs/libBATmvc.rootmap: libBATmvc.la
	cp libBATmvc.rootmap libBATmvc_rdict.pcm .libs/

.libs/libBATmvc_rdict.pcm: .libs/libBATmvc.rootmap

else # not WITH_CLING
### ROOT-5: ###

libBATmvc_rdict.cxx: $(library_include_HEADERS) LinkDef.h
	$(ROOTCINT) -f $@ -c -p $(CPPFLAGS) $(CXXFLAGS) -I.includes -I$(top_builddir)/src/.includes $+

.libs/libBATmvc.rootmap: libBATmvc.la | $(library_include_HEADERS) LinkDef.h
	$(RLIBMAP) -o $@.tmp -l .libs/libBATmvc@SHLIBEXT@ -c $|
	$(SED) 's|\.so$$|'@SHLIBEXT@'|' $@.tmp > $@ && $(RM) $@.tmp

endif # WITH_CLING
