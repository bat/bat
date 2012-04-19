TEXFILE=quickstart
EPSFILES=bat.eps

.PHONY : latex html mostly.clean clean

all: latex html

latex: $(TEXFILE).pdf
latex: $(TEXFILE).html

%.dvi %.aux: %.tex
	latex $< && latex $< && latex $<

%.ps: %.dvi $(EPSFILES)
	dvips $< -o $@

%.pdf : %.ps
	dvipdf -sPAPERSIZE=a4 -dPDFSETTINGS=/prepress $(*F)

%.html : %.tex %.aux
	./tth $<

mostly.clean:
	rm -f $(addprefix $(TEXFILE),.toc .dvi .aux .log)

clean: mostly.clean
	rm -f $(addprefix $(TEXFILE),.pdf .ps .html)