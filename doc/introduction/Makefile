TEXFILE = introduction

.PHONY : latex clean

all: latex

latex: $(TEXFILE).pdf

%.pdf: %.tex
	pdflatex $<
	pdflatex $<
	pdflatex $<

clean:
	rm -f $(addprefix $(TEXFILE), .aux .dvi .log .pdf .out .toc)
