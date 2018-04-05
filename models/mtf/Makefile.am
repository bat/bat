CLEANFILES = *_rdict.* *_rdict_* *.rootmap

clean-local:
	-rm -rf .includes

lib_LTLIBRARIES= libBATmtf.la
libBATmtf_la_LDFLAGS= -version-info 1:0:0 $(ROOTLIBDIRFLAGS)
libBATmtf_la_LIBADD = $(top_builddir)/src/libBAT.la $(ROOTLIBS)

AM_CXXFLAGS = $(ROOTAUXCFLAGS) -Wall -fPIC -I.includes -I$(top_builddir)/src/.includes $(ROOTCFLAGS)
AM_CXXFLAGS += @AM_CXXFLAGS@
AM_LDFLAGS  = -O
AM_SOFLAGS  = -shared

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

EXTRA_DIST = $(incs)

nodist_library_include_HEADERS = $(incs:%.h=.includes/$(PACKAGE)/%.h)

$(nodist_library_include_HEADERS): .includes/$(PACKAGE)/%.h: %.h
	mkdir -p .includes/$(PACKAGE)
	sed 's|include \"\(BC.*.h\)\"|include <'$(PACKAGE)'\/\1>| ; s/include \"\.\.\/\.\.\/\('$(PACKAGE)'\/BC.*.h\)\"/include <\1>/' $< > $@

libBATmtf_la_SOURCES= \
	$(incs:.h=.cxx)

nodist_libBATmtf_la_SOURCES= \
	libBATmtf_rdict.cxx

BUILT_SOURCES = libBATmtf_rdict.cxx

EXTRA_DIST += LinkDef.h

install-exec-hook:
	$(CREATE_SO_LINK)

rootmapdir = $(libdir)
rootmap_DATA = .libs/libBATmtf.rootmap

if WITH_CLING
### ROOT-6: ###

rootmap_DATA += .libs/libBATmtf_rdict.pcm

# short hand to make rule more concise, rtmap is not a standard automake variable
rtmap = libBATmtf.rootmap
libBATmtf_rdict.cxx: $(nodist_library_include_HEADERS) LinkDef.h
	$(ROOTCLING) -f $@.tmp -s libBATmtf@SHLIBEXT@ -rml libBATmtf@SHLIBEXT@ -rmf libBATmtf.rootmap -c $(CPPFLAGS) $(CXXFLAGS) -I.includes -I$(top_builddir)/src/.includes -I$(includedir) $+
	@# Some magic to prefix header names with "$(PACKAGE)/", and only that, in dictionary and rootmap:
	@# avoid sed -i: not portable between Mac OS X and linux => use back-up files
	$(GREP) -F -v '"'"`pwd`"'/",' $@.tmp | grep -v '^\"\([^/"]*/\)*\.includes\",*$$' | $(SED) 's|"\([^/"]*/\)*\([^/"]*[.]h\)",|"'$(PACKAGE)/'\2",| ; s|\\"\([^/"]*/\)*\([^/"]*[.]h\)\\"\\n"|\\"'$(PACKAGE)/'\2\\"\\n"|' > $@.tmp2
	$(SED) 's|\$$clingAutoload\$$\([^/"]*/\)*|$$clingAutoload$$'$(PACKAGE)'/|; /.*DICTPAYLOAD(.*/,/.*)DICTPAYLOAD.*/ s|#include "\([^/"]*/\)*\(.*\)"|#include <'$(PACKAGE)'/\2>|' $@.tmp2 > $@ && $(RM) $@.tmp $@.tmp2
	$(SED) 's|\(header  *\)\([^ ].*/\)*\([^ ].*[.]h\)|\1'$(PACKAGE)/'\3|' $(rtmap) > $(rtmap).tmp && mv $(rtmap).tmp $(rtmap)

.libs/libBATmtf.rootmap: libBATmtf.la
	cp libBATmtf.rootmap libBATmtf_rdict.pcm .libs/

.libs/libBATmtf_rdict.pcm: .libs/libBATmtf.rootmap

else # not WITH_CLING
### ROOT-5: ###

libBATmtf_rdict.cxx: $(nodist_library_include_HEADERS) LinkDef.h
	$(ROOTCINT) -f $@ -c -p $(CPPFLAGS) $(CXXFLAGS) -I.includes -I$(top_builddir)/src/.includes $+

.libs/libBATmtf.rootmap: libBATmtf.la | $(nodist_library_include_HEADERS) LinkDef.h
	$(RLIBMAP) -o $@.tmp -l .libs/libBATmtf@SHLIBEXT@ -c $|
	$(SED) 's|\.so$$|'@SHLIBEXT@'|' $@.tmp > $@ && $(RM) $@.tmp

endif # WITH_CLING
