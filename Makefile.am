AUTOMAKE_OPTIONS = foreign
ACLOCAL = aclocal
AUTOMAKE= automake
clean-local:
	rm -rf autom4te.cache
SUBDIRS= \
	src
EXTRA_DIST= \
	src/libBAT.rootmap \
	doc/COPYING doc/LICENSE doc/changelog doc/releasenotes.txt \
	examples
dist-hook: prep-docs prep-examples
	rm -f $(distdir)/src/*Dict.*
	rm -rf `find $(distdir)/examples -name CVS`
prep-docs: BAT-intro.ps
BAT-intro.ps:
	mkdir intro-tmp
	cp -r doc/introduction intro-tmp
	cd intro-tmp/introduction; latex introduction; latex introduction; dvips introduction
	mv intro-tmp/introduction/introduction.ps $(distdir)/doc/BAT-intro.ps
	rm -rf intro-tmp


