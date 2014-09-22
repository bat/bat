TEXFILE = introduction
EPSFILES = bat.eps

.PHONY : latex clean

all: latex

latex: $(TEXFILE).pdf

%.pdf: %.tex
	pdflatex $<
	pdflatex $<

clean:
	rm -f $(addprefix $(TEXFILE), .aux .dvi .log .pdf .out .toc)
