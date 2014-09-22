AUTOMAKE_OPTIONS = foreign
ACLOCAL = aclocal
AUTOMAKE= automake

clean-local:
	rm -rf autom4te.cache

SUBDIRS= \
	src \
	models/base \
	models/mtf \
	models/mvc

EXTRA_DIST= \
	INSTALL.md \
	README.md \
	src/libBAT.rootmap \
	models/base/BCRooInterface.rootmap \
	models/base/libBATmodels.rootmap.in \
	models/base/BATCalculatorLinkDef.h \
	models/mtf/libBATmtf.rootmap \
	models/mvc/libBATmvc.rootmap \
	doc/COPYING doc/CREDITS doc/LICENSE doc/ChangeLog doc/releasenotes.txt \
	examples \
	tools/CreateProject.sh.in

dist_bin_SCRIPTS = bat-config

dist-hook: prep-docs
#	rm -f $(distdir)/src/*Dict.*
#	rm -f $(distdir)/models/*Dict.*
	rm -f $(distdir)/models/libBATmodels.rootmap
	rm -rf `find $(distdir)/examples -name Makefile`

prep-docs: BAT-intro.pdf BAT-quickstart.pdf

BAT-intro.pdf:
	cd $(top_srcdir)/doc/introduction \
	make \
	mv introduction.pdf $(distdir)/doc/BAT-intro.pdf

BAT-quickstart.pdf:
	cd $(top_srcdir)/doc/quickstart \
	make \
	mv quickstart.pdf $(distdir)/doc/BAT-quickstart.pdf

if COND_PKGCONFIG
pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = bat.pc
endif
