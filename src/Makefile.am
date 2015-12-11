CLEANFILES = *_rdict.* *_rdict_* *.rootmap

clean-local:
	-rm -rf .includes

lib_LTLIBRARIES = libBAT.la
libBAT_la_LDFLAGS = -version-info 1:0:0

# this fixes the linker error on ubuntu with gcc > 4.4 with missing reference to TCanvas(bool) etc.
libBAT_la_LDFLAGS += -L$(ROOTLIBDIR)
libBAT_la_LIBADD = $(ROOTLIBS)

AM_CXXFLAGS = $(ROOTAUXCFLAGS) -Wall -fPIC -I$(top_srcdir)/BAT -I.includes $(ROOTCFLAGS) $(CUBAINCLUDE)
AM_LDFLAGS  = -O
AM_SOFLAGS  = -shared

if CUBA
libBAT_la_LIBADD += $(CUBALIB)
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
	BCHistogramBase.h \
	BCH1D.h \
	BCH2D.h \
	BCPrior.h \
	BCConstantPrior.h \
	BCGaussianPrior.h \
	BCSplitGaussianPrior.h \
	BCCauchyPrior.h \
	BCTF1Prior.h \
	BCTF1LogPrior.h \
	BCTH1Prior.h \
	BCPositiveDefinitePrior.h \
	BCObservableSet.h \
	BCParameterSet.h \
	BCObservable.h \
	BCParameter.h \
	BCVariable.h \
	BCEngineMCMC.h \
	BCIntegrate.h \
	BCModel.h \
	BCEmptyModel.h \
	BCPriorModel.h \
	BCModelManager.h \
	BCLog.h \
	BCMath.h \
	BCAux.h

non_source_incs = \
	BCVariableSet.h

EXTRA_DIST += $(incs:%.h=$(top_srcdir)/BAT/%.h) $(non_source_incs:%.h=$(top_srcdir)/BAT/%.h)

nodist_library_include_HEADERS = $(incs:%.h=.includes/$(PACKAGE)/%.h) $(non_source_incs:%.h=.includes/$(PACKAGE)/%.h)

$(nodist_library_include_HEADERS): .includes/$(PACKAGE)/%.h: $(top_srcdir)/BAT/%.h
	mkdir -p .includes/$(PACKAGE)
	sed 's|include \"\(BC.*.h\)\"|include <'$(PACKAGE)'\/\1>| ; s/include \"\.\.\/\.\.\/\('$(PACKAGE)'\/BC.*.h\)\"/include <\1>/' $< > $@

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
		$(SED) -e 's/include \"\(BC.*.h\)\"/include <BAT\/\1>/' $$i > src-tmp/`basename $$i`.h; \
		if ! diff -q src-tmp/`basename $$i`.h "$$i" > /dev/null; then \
			echo "updating $$i"; \
			cp -f src-tmp/`basename $$i`.h "$$i"; \
			chmod 664 "$$i"; \
		fi; \
	done; \
	rm -rf src-tmp;

# rootmapdir implicitly used in rootmap_DATA
rootmapdir = $(libdir)
rootmap_DATA = .libs/libBAT.rootmap

if WITH_CLING
### ROOT-6: ###

rootmap_DATA += .libs/libBAT_rdict.pcm

# short hand to make rule more concise, rtmap is not a standard automake variable
rtmap = libBAT.rootmap
libBAT_rdict.cxx: $(nodist_library_include_HEADERS) LinkDef.h
	$(ROOTCLING) -f $@.tmp -s libBAT@SHLIBEXT@ -rml libBAT@SHLIBEXT@ -rmf libBAT.rootmap -c $(CPPFLAGS) $(CXXFLAGS) -I.includes -I$(includedir) $+
	@# Some magic to prefix header names with "$(PACKAGE)/", and only that, in dictionary and rootmap:
	@# avoid sed -i: not portable between Mac OS X and linux => use back-up files
	$(GREP) -F -v '"'"`pwd`"'/",' $@.tmp | grep -v '^\"\([^/"]*/\)*\.includes\",*$$' | $(SED) 's|"\([^/"]*/\)*\([^/"]*[.]h\)",|"'$(PACKAGE)/'\2",| ; s|\\"\([^/"]*/\)*\([^/"]*[.]h\)\\"\\n"|\\"'$(PACKAGE)/'\2\\"\\n"|' > $@.tmp2
	$(SED) 's|\$$clingAutoload\$$\([^/"]*/\)*|$$clingAutoload$$'$(PACKAGE)'/|; /.*DICTPAYLOAD(.*/,/.*)DICTPAYLOAD.*/ s|#include "\([^/"]*/\)*\(.*\)"|#include <'$(PACKAGE)'/\2>|' $@.tmp2 > $@ && $(RM) $@.tmp $@.tmp2
	$(SED) 's|\(header  *\)\([^ ].*/\)*\([^ ].*[.]h\)|\1'$(PACKAGE)/'\3|' $(rtmap) > $(rtmap).tmp && mv $(rtmap).tmp $(rtmap)

.libs/libBAT.rootmap: libBAT.la
	cp libBAT.rootmap libBAT_rdict.pcm .libs/

.libs/libBAT_rdict.pcm: .libs/libBAT.rootmap

else # not WITH_CLING
### ROOT-5: ###

libBAT_rdict.cxx: $(nodist_library_include_HEADERS) LinkDef.h
	$(ROOTCINT) -f $@ -c -p $(CPPFLAGS) $(CXXFLAGS) -I.includes $+

.libs/libBAT.rootmap: libBAT.la | $(nodist_library_include_HEADERS) LinkDef.h
	$(RLIBMAP) -o $@.tmp -l .libs/libBAT@SHLIBEXT@ -c $|
	$(SED) 's|\.so$$|'@SHLIBEXT@'|' $@.tmp > $@ && $(RM) $@.tmp

endif # WITH_CLING
