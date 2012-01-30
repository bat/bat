AUTOMAKE_OPTIONS = foreign
ACLOCAL = aclocal
AUTOMAKE= automake

clean-local:
	rm -rf autom4te.cache

SUBDIRS= \
	src \
	models \
	benchmarks/test

EXTRA_DIST= \
	src/libBAT.rootmap \
	models/BCRooInterface.rootmap \
	models/libBATmodels.rootmap.in \
	models/BATCalculatorLinkDef.h \
	doc/COPYING doc/CREDITS doc/LICENSE doc/ChangeLog doc/releasenotes.txt \
	examples \
	tools/CreateProject.sh.in

dist-hook: prep-docs
#	rm -f $(distdir)/src/*Dict.*
#	rm -f $(distdir)/models/*Dict.*
	rm -f $(distdir)/models/libBATmodels.rootmap
	rm -rf `find $(distdir)/examples -name .svn`
	rm -rf `find $(distdir)/models -name .svn`
	rm -rf `find $(distdir)/examples -name Makefile`

prep-docs: BAT-intro.ps

BAT-intro.ps:
	mkdir intro-tmp
	cp -r doc/introduction intro-tmp
	cd intro-tmp/introduction; latex introduction; latex introduction; dvips -o introduction.ps introduction
	mv intro-tmp/introduction/introduction.ps $(distdir)/doc/BAT-intro.ps
	rm -rf intro-tmp


