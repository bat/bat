# every file is a chapter of the manual; the order of files is respected
input = \
		getting-started.md \
		bayes.md \
		integration.md \
		mcmc.md \
		optimization.md \
		predefined-models.md \
		mtf.md \
		\
		code-structure.md \
		model-comparison.md \
		output.md \
		factorized-prior.md \
		sharing-samples.md \
		threading.md
# Chapters temporarily deactivated because they are empty
#   validation.md

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
	IMAGE_PATH = $(srcdir)/.. $(srcdir)/figures $(srcdir)/examples $(top_builddir)/examples/basic/tutorial
	EXAMPLE_PATH = $(srcdir)/examples $(top_srcdir)/BAT
	CITE_BIB_FILES = $(srcdir)/references.bib
	EXTRA_PACKAGES = bat
endef
export extradoxyconf

# do nothing by default so `make` in top-level dir doesn't require all the
# tools used for the docs. To build both html and pdf manual, use the `pdf` target
all: ;

####
# 1. run tutorial example to create pdf plots for latex manual
# 2. postprocess: extract interesting page out of multipage pdf
# 3. then convert to png for html manual
####

tut_dir = $(top_builddir)/examples/basic/tutorial
tut_file = $(tut_dir)/gaus_mod_

# Let pdfs depend on log file to capture the build dependency that the example
# creates multiple plots. This has a weakness in incremental builds if the log
# file is present but one of the pdfs isn't. But it's safe in a from-scratch
# build. See
# https://www.cmcrossroads.com/article/rules-multiple-outputs-gnu-make
$(tut_dir)/test-suite.log:
	$(MAKE) -C $(top_builddir)/examples/basic/tutorial check

tutorial_plots= \
	$(tut_file)correlation.pdf \
	$(tut_file)correlationMatrix.pdf \
	$(tut_file)parameters.pdf \
	$(tut_file)plots.pdf \
	$(tut_file)update.pdf

$(tutorial_plots): $(tut_dir)/test-suite.log

# Reduce typing here by defining a function
# extract_page PAGE INPUT OUTPUT
define extract_page
	pdftk $(2) cat $(1) output $(3)
endef

$(tut_file)parameters-1.pdf: $(tut_file)parameters.pdf
	$(call extract_page, 1, $?, $@)
$(tut_file)plots-1.pdf: $(tut_file)plots.pdf
	$(call extract_page, 1, $?, $@)
$(tut_file)plots-2.pdf: $(tut_file)plots.pdf
	$(call extract_page, 3, $?, $@)
$(tut_file)update-1.pdf: $(tut_file)update.pdf
	$(call extract_page, 1, $?, $@)

tutorial_plots_png = \
	$(tut_file)correlation.png \
	$(tut_file)correlationMatrix.png \
	$(tut_file)parameters-1.png \
	$(tut_file)plots-1.png \
	$(tut_file)plots-2.png \
	$(tut_file)update-1.png

DOXYDEPENDENCIES = $(tutorial_plots_png)

# convert pdf to png and crop
%.png: %.pdf
	pdftoppm $? $(basename $@) -png -r 100 -cropbox -singlefile

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
