input = front.md basics.md bayes.md mcmc.md factorized.md share.md
doxyinput = $(addprefix $(srcdir)/, $(input))
stylefile = bat.sty

define extradoxyconf
	TAGFILES = $(top_builddir)/doc/ref-guide/bat-ref.tag=$(top_builddir)/ref-guide/html
	IMAGE_PATH = $(srcdir)/figures
	CITE_BIB_FILES = $(srcdir)/references.bib
endef
export extradoxyconf

include $(top_srcdir)/doc/Makefile-common.am

EXTRA_DIST = \
	$(input) Doxyfile references.bib $(stylefile) \
	figures/random-walk.png figures/random-walk.pdf

# rely on doxygen to create latex/
.PHONY: pdf
pdf: doxy
	cp -f $(srcdir)/$(stylefile) $(builddir)/latex/
	make -C latex $@

.PHONY: print
print:
	$(info $(doxyconf))
	$(info $(extradoxyconf))
	@echo
	@echo "$$doxyconf"
	@echo "$$extradoxyconf"
