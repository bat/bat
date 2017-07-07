Building the documentation
==========================

To build only reference guide with `doxygen`, do

    make ref-guide

then open `ref-guide/html/index.html` with a web browser.

To build the manual, the reference guide is built first as a
dependency when calling

    make manual

The pdf and html outputs are in `manual/BAT-manual.pdf` and in
`manual/html/index.html`. To select only one output format, do

    make -C manual doxy # html
    make -C manual pdf

The manual makes use of math expressions rendered with mathjax in the
html output. To increase the rendering speed, one can use a local copy
whose path needs to be defined before the manual is built.

    export MATHJAX_RELPATH=/usr/share/javascript/mathjax/
    make manual

Extra dependencies
------------------

One needs `doxygen, graphviz` and `latex`. For latex, the metapackage
`texlive-collection-latexrecommended` on cent os has what is needed.
