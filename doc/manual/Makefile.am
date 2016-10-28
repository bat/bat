input = $(addprefix $(top_srcdir)/doc/manual/, front.md basics.md bayes.md mcmc.md factorized.md share.md)

define extradoxyconf
	TAGFILES = $(top_builddir)/doc/ref-guide/bat-ref.tag=$(top_builddir)/ref-guide/html
endef
export extradoxyconf

include $(top_srcdir)/doc/Makefile-common.am

.PHONY: pdf
pdf: doxy
	cp -f $(top_srcdir)/doc/manual/bat.sty $(top_builddir)/doc/manual/latex
	make -C latex $@

.PHONY: print
print:
	$(info $(doxyconf))
	$(info $(extradoxyconf))
	@echo
	@echo "$$doxyconf"
	@echo "$$extradoxyconf"
