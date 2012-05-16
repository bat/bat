AUTOMAKE_OPTIONS = foreign
ACLOCAL = aclocal
AUTOMAKE= automake

clean-local:
	rm -rf autom4te.cache

SUBDIRS= \
	src \
	models/base \
	models/mtf

EXTRA_DIST= \
	src/libBAT.rootmap \
	models/base/BCRooInterface.rootmap \
	models/base/libBATmodels.rootmap.in \
	models/base/BATCalculatorLinkDef.h \
	models/mtf/libBATmtf.rootmap \
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

prep-docs: BAT-intro.pdf BAT-quickstart.pdf

BAT-intro.pdf:
	mkdir intro-tmp
	cp -r doc/introduction intro-tmp
	cd intro-tmp/introduction; latex introduction; latex introduction; dvips -o introduction.ps introduction; ps2pdf introduction.ps
	mv intro-tmp/introduction/introduction.pdf $(distdir)/doc/BAT-intro.pdf
	rm -rf intro-tmp

BAT-quickstart.pdf:
	mkdir quickstart-tmp
	cp -r doc/quickstart quickstart-tmp
	cd quickstart-tmp/quickstart; latex quickstart; latex quickstart; dvips -o quickstart.ps quickstart; ps2pdf quickstart.ps
	mv quickstart-tmp/quickstart/quickstart.pdf $(distdir)/doc/BAT-quickstart.pdf
	rm -rf quickstart-tmp

