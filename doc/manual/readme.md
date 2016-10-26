create the reference guide first to create a tag file that is referenced from the manual

To make latex recognize `\newcommand`, put a symbolic link to `bat.sty` into `latex` subdirectory, then call `make pdf` there.

To define the same commands also for mathjax, follow http://stackoverflow.com/questions/40270302
