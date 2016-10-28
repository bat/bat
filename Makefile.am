AUTOMAKE_OPTIONS = foreign
ACLOCAL  = aclocal
AUTOMAKE = automake

clean-local:
	-rm -rf autom4te.cache
	-rm -rf external/cuba-*

# files build by configure on the user side
DISTCLEANFILES = \
	bat.pc \
	tools/bat-project \
	doc/introduction/versions.tex

SUBDIRS= \
	src \
	models/base \
	models/mtf \
	test \
	doc \
	examples

# TODO remove intro, quickstart
EXTRA_DIST= \
	INSTALL.md \
	README.md \
	models/base/LDef.h  \
	models/base/LDef_roostats.h \
	doc/introduction/bat.pdf  doc/introduction/introduction.tex  doc/introduction/Makefile \
	doc/quickstart/bat.pdf  doc/quickstart/Makefile  doc/quickstart/quickstart.tex \
	examples \
	tools/bat-project.in

dist_bin_SCRIPTS = \
	bat-config \
	tools/bat-project

# Must not fix the installation dir at configure time,
# so we have to create the file. For the reasoning, see
# https://www.gnu.org/software/autoconf/manual/autoconf-2.63/html_node/Installation-Directory-Variables.html#Installation-Directory-Variables
edit = sed \
	-e 's|@datadir[@]|$(pkgdatadir)|g' \
	-e 's|@prefix[@]|$(prefix)|g'

tools/bat-project: Makefile $(srcdir)/tools/bat-project.in
	rm -f $@ $@.tmp
	$(edit) '$(srcdir)/$@.in' >$@.tmp
	chmod +x $@.tmp
	chmod a-w $@.tmp
	mv $@.tmp $@

dist_pkgdata_DATA = \
	tools/Main1Template.cxx \
	tools/Main2Template.cxx \
	tools/ModelTemplate.h tools/ModelTemplate.cxx \
	tools/MakefileTemplate

dist-hook: prep-docs
	rm -f $(distdir)/models/libBATmodels.rootmap
	rm -rf `find $(distdir)/examples -name Makefile`

# TODO remove or move to doc/Makefile.am
prep-docs: BAT-intro.pdf BAT-quickstart.pdf

BAT-intro.pdf:
	mkdir -p $(top_builddir)/doc/intro-tmp
	cp -r $(top_srcdir)/doc/introduction/* $(top_builddir)/doc/introduction/versions.tex $(top_builddir)/doc/intro-tmp
	cd $(top_builddir)/doc/intro-tmp && $(MAKE) clean && $(MAKE)
	mv $(top_builddir)/doc/intro-tmp/introduction.pdf $(distdir)/doc/BAT-intro.pdf
	rm -rf $(top_builddir)/doc/intro-tmp

BAT-quickstart.pdf:
	mkdir -p $(top_builddir)/doc/quick-tmp
	cp -r $(top_srcdir)/doc/quickstart/* $(top_builddir)/doc/quick-tmp
	cd $(top_builddir)/doc/quick-tmp && $(MAKE) clean && $(MAKE)
	mv $(top_builddir)/doc/quick-tmp/quickstart.pdf $(distdir)/doc/BAT-quickstart.pdf
	rm -rf $(top_builddir)/doc/quick-tmp

if COND_PKGCONFIG
pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = bat.pc
endif

if COND_USE_RPATH
USE_RPATH=yes
else
USE_RPATH=no
endif


install-data-hook:
	@echo
	@echo "   ============================================================"
	@echo "                   BAT INSTALLATION SUCCESSFUL!          "
	@echo
	@echo "In the bash, you can set up your environment with"
	@echo
	@echo 'export PATH="$(bindir):$$PATH"'
	@test "$(USE_RPATH)" != "yes" && echo 'export LD_LIBRARY_PATH="$(libdir):$$LD_LIBRARY_PATH"' || true
	@echo 'export CPATH="$(includedir):$$CPATH"'
	@echo 'export PKG_CONFIG_PATH="$(pkgconfigdir):$$PKG_CONFIG_PATH"'
	@echo
	@echo "   ============================================================"
	@echo

# format all C++ source files. Exclude certain files for which astyle
# doesn't work.
.PHONY: format
format:
	cd $(top_srcdir) && \
	astyle --options=tools/astylerc --formatted `find BAT/ benchmarks/ examples/ tools $(SUBDIRS) -name '*.C' -o -name '*.cxx' -o -name '*.h' | egrep -v -f tools/astyle-ignore`
