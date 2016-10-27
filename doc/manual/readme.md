Build
-----
create the reference guide first to create a tag file that is referenced from the manual

To make latex recognize `\newcommand`, put a symbolic link to `bat.sty` into `latex` subdirectory, then call `make pdf` there.

To define the same commands also for mathjax, follow http://stackoverflow.com/questions/40270302

For deployment, change

    MATHJAX_RELPATH        = http://cdn.mathjax.org/mathjax/latest

from a local version

adding a chapter
----------------

Create file 'chapter3.md', add in the proper location in `doxyfile`

     INPUT = front.md basics.md bayes.md chapter3.md

Copy over the labeling of sections from `bayes.md`, it doesn't match the doxygen docs on markdown. Add `[TOC]`!

debugging doxygen
-----------------

    doxygen Doxyfile -d

images
------

From inkscape, export to png for html output with

    inkscape --export-area-drawing --export-png=random-walk.png --export-width 500 random-walk.svg

The markdown image support is not good enough for us, so import images separately for html and latex. Forcibly add both images to git!

Keep the caption in `""` short because it doesn't look good in HTML and avoid math as mathjax fails in caption. Unfortunately caption has to be copied

    @anchor random-walk-2D
    @image html random-walk.png "2D random walk."
    @image latex random-walk.pdf "2D random walk." width=0.5\textwidth

Refer to the image with

    @ref random-walk-2D "2D example plot"

output specific code
----------

    @htmlonly
    some text only visible in the html output
    @endhtmlonly

    @latexonly
    some text only visible in the latex output
    @endlatexonly
