CLEANFILES = *_rdict.* *_rdict_* *.rootmap LinkDef.h

clean-local:
	-rm -rf .includes

lib_LTLIBRARIES= libBATmodels.la
libBATmodels_la_LDFLAGS= -version-info 1:0:0 $(ROOTLIBDIRFLAGS)
libBATmodels_la_LIBADD = $(top_builddir)/src/libBAT.la $(ROOTLIBS)

AM_CXXFLAGS = $(ROOTAUXCFLAGS) -Wall -fPIC -I.includes -I$(top_builddir)/src/.includes $(ROOTCFLAGS)
AM_CXXFLAGS += @AM_CXXFLAGS@
AM_LDFLAGS  = -O
AM_SOFLAGS  = -shared

if IS_DARWIN
CREATE_SO_LINK = cd $(DESTDIR)$(libdir); rm -f libBATmodels.so; ln -s libBATmodels.dylib libBATmodels.so
endif

library_includedir= $(includedir)/$(PACKAGE)

incs = \
	BCFitter.h \
	BCEfficiencyFitter.h \
	BCGraphFitter.h \
	BCHistogramFitter.h

roostat_incs = \
	BCRooInterface.h \
	BATCalculator.h

linkdefs = LDef.h

## RooStats specific options
if ROOSTATS
incs += $(roostat_incs)
linkdefs += LDef_roostats.h
endif
## End of RooStats specific options

LinkDef.h: $(linkdefs)
	cat $+ > $@

EXTRA_DIST = $(incs) $(roostat_incs)

nodist_library_include_HEADERS = $(incs:%.h=.includes/$(PACKAGE)/%.h)

$(nodist_library_include_HEADERS): .includes/$(PACKAGE)/%.h: %.h
	mkdir -p .includes/$(PACKAGE)
	sed 's|include \"\(BC.*.h\)\"|include <'$(PACKAGE)'\/\1>| ; s/include \"\.\.\/\.\.\/\('$(PACKAGE)'\/BC.*.h\)\"/include <\1>/' $< > $@

libBATmodels_la_SOURCES= \
	$(incs:.h=.cxx)

nodist_libBATmodels_la_SOURCES= \
	libBATmodels_rdict.cxx

BUILT_SOURCES = libBATmodels_rdict.cxx LinkDef.h

EXTRA_DIST += LDef.h LDef_roostats.h

install-exec-hook:
	$(CREATE_SO_LINK)

rootmapdir = $(libdir)
rootmap_DATA = .libs/libBATmodels.rootmap

if WITH_CLING
### ROOT-6: ###

rootmap_DATA += .libs/libBATmodels_rdict.pcm

# short hand to make rule more concise, map is not a standard automake variable
rtmap = libBATmodels.rootmap
libBATmodels_rdict.cxx: $(nodist_library_include_HEADERS) LinkDef.h
	$(ROOTCLING) -f $@.tmp -s libBATmodels@SHLIBEXT@ -rml libBATmodels@SHLIBEXT@ -rmf libBATmodels.rootmap -c $(CPPFLAGS) $(CXXFLAGS) -I.includes -I$(top_builddir)/src/.includes -I$(includedir) $+
	@# Some magic to prefix header names with "$(PACKAGE)/", and only that, in dictionary and rootmap
	@# avoid sed -i: not portable between Mac OS X and linux => use back-up files
	$(GREP) -F -v '"'"`pwd`"'/",' $@.tmp | grep -v '^\"\([^/"]*/\)*\.includes\",*$$' | $(SED) 's|"\([^/"]*/\)*\([^/"]*[.]h\)",|"'$(PACKAGE)/'\2",| ; s|\\"\([^/"]*/\)*\([^/"]*[.]h\)\\"\\n"|\\"'$(PACKAGE)/'\2\\"\\n"|' > $@.tmp2
	$(SED) 's|\$$clingAutoload\$$\([^/"]*/\)*|$$clingAutoload$$'$(PACKAGE)'/|; /.*DICTPAYLOAD(.*/,/.*)DICTPAYLOAD.*/ s|#include "\([^/"]*/\)*\(.*\)"|#include <'$(PACKAGE)'/\2>|' $@.tmp2 > $@ && $(RM) $@.tmp $@.tmp2
	$(SED) 's|\(header  *\)\([^ ].*/\)*\([^ ].*[.]h\)|\1'$(PACKAGE)/'\3|' $(rtmap) > $(rtmap).tmp && mv $(rtmap).tmp $(rtmap)

.libs/libBATmodels.rootmap: libBATmodels.la
	cp libBATmodels.rootmap libBATmodels_rdict.pcm .libs/

.libs/libBATmodels_rdict.pcm: .libs/libBATmodels.rootmap

else # not WITH_CLING
### ROOT-5: ###

libBATmodels_rdict.cxx: $(nodist_library_include_HEADERS) LinkDef.h
	$(ROOTCINT) -f $@ -c -p $(CPPFLAGS) $(CXXFLAGS) -I.includes -I$(top_builddir)/src/.includes $+

.libs/libBATmodels.rootmap: libBATmodels.la | $(nodist_library_include_HEADERS) LinkDef.h
	$(RLIBMAP) -o $@.tmp -l .libs/libBATmodels@SHLIBEXT@ -c $|
	$(SED) 's|\.so$$|'@SHLIBEXT@'|' $@.tmp > $@ && $(RM) $@.tmp

endif # WITH_CLING
