default: latex

latex: 
	latex introduction.tex
	latex introduction.tex
	dvips introduction.dvi -o introduction.ps
	dvipdf -sPAPERSIZE=a4 -dPDFSETTINGS=/prepress introduction

clean:
	rm -f *.toc *.dvi *.aux