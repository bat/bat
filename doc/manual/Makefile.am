# every file is a chapter of the manual; the order of files is respected
input = \
		getting-started.md \
		bayes.md \
		mcmc.md \
		integration.md \
		optimization.md \
		predefined-models.md \
		\
		code-structure.md \
		output.md \
		model-comparison.md \
		factorized-prior.md \
		sharing-samples.md \
		threading.md \
		validation.md \
		troubleshooting.md

doxyinput = $(srcdir)/front.md
doxyinput += $(top_srcdir)/INSTALL.md
doxyinput += $(input:%=$(srcdir)/%)
stylefile = bat.sty

define extradoxyconf
	PROJECT_NAME = "BAT manual"
	USE_MDFILE_AS_MAINPAGE = front.md
	LAYOUT_FILE = $(srcdir)/DoxygenLayout.xml
	TAGFILES = $(top_builddir)/doc/ref-guide/bat-ref.tag=$(top_builddir)/ref-guide/html
	GENERATE_LATEX = YES
	IMAGE_PATH = $(srcdir)/figures
	CITE_BIB_FILES = $(srcdir)/references.bib
	EXTRA_PACKAGES = bat
endef
export extradoxyconf

include $(top_srcdir)/doc/Makefile-common.am

EXTRA_DIST = \
	$(input) Doxyfile DoxygenLayout.xml \
	references.bib $(stylefile) \
	figures/random-walk.png figures/random-walk.pdf

# rely on doxygen to create latex/
.PHONY: pdf
pdf: doxy
	cp -f $(srcdir)/$(stylefile) $(builddir)/latex/
	make -C latex $@
	mv latex/refman.pdf @BAT_MANUAL_PDF@

clean-local:
	-rm -rf doxygen_sqlite3.db html/ latex/ @BAT_MANUAL_PDF@
