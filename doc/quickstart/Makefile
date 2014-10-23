TEXFILE = quickstart

.PHONY : latex clean

all: latex

latex: $(TEXFILE).pdf

%.pdf: %.tex
	pdflatex $<
	pdflatex $<

clean:
	rm -f $(addprefix $(TEXFILE), .aux .dvi .log .out .pdf .toc)
