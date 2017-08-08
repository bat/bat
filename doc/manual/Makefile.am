# every file is a chapter of the manual; the order of files is respected
input = \
		getting-started.md \
		bayes.md \
		integration.md \
		mcmc.md \
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

TAGFILE = bat-ref.tag

define extradoxyconf
	PROJECT_NAME = "BAT manual"
	USE_MDFILE_AS_MAINPAGE = front.md
	LAYOUT_FILE = $(srcdir)/DoxygenLayout.xml
	TAGFILES = $(top_builddir)/doc/ref-guide/$(TAGFILE)=$(top_builddir)/ref-guide/html
	GENERATE_LATEX = YES
	IMAGE_PATH = $(srcdir)/figures $(srcdir)/examples
	EXAMPLE_PATH = $(srcdir)/examples $(top_srcdir)/BAT
	CITE_BIB_FILES = $(srcdir)/references.bib
	EXTRA_PACKAGES = bat
endef
export extradoxyconf

all: pdf

# for doxygen to succeed the TAGFILE has to exist
include $(top_srcdir)/doc/Makefile-common.am

EXTRA_DIST = \
	$(input) DoxygenLayout.xml \
	references.bib $(stylefile) \
	figures/random-walk.png figures/random-walk.pdf \
	figures/histogram.png figures/histogram.pdf

.PHONY: $(TAGFILE)
$(TAGFILE):
	$(MAKE) -C ../ref-guide $@

# doxygen creates latex sources but doesn't build pdf.
.PHONY: pdf
pdf: doxy
	cp -f $(srcdir)/$(stylefile) $(builddir)/latex/
	$(MAKE) -C latex $@
	mv latex/refman.pdf @BAT_MANUAL_PDF@

clean-local:
	-rm -rf doxygen_sqlite3.db html/ latex/ @BAT_MANUAL_PDF@
