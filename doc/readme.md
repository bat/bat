Building the documentation
==========================

To build the entire documentation, do

    make

To build only the reference guide with `doxygen`, do

    make ref-guide doxy

then open `ref-guide/html/index.html` with a web browser.

To build the manual in html and pdf format, the reference guide is built first
as a dependency when calling

    make -C manual pdf
    make -C manual doxy

The pdf and html outputs are in `manual/BAT-manual.pdf` and in
`manual/html/index.html`.

The manual makes use of math expressions rendered with mathjax in the
html output. To increase the rendering speed, one can use a local copy
whose path needs to be defined before the manual is built.

    export MATHJAX_RELPATH=/usr/share/javascript/mathjax/
    make -C manual

Extra dependencies
------------------

One needs `doxygen, graphviz, pdftk, pdftoppm` and `pdflatex`. For latex, the metapackage
`texlive-collection-latexrecommended` on cent os has what is needed.
