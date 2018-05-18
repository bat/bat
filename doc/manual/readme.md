Build
-----


### high level

follow the instructions at `doc/readme.md`:

### details

The reference guide needs to first be created and it produces a tag file that is referenced from the manual. The makefile handles this dependency automatically.

To make latex recognize `\newcommand`, add it to `doc/bat.sty`, it is copied by the makefile to the latex subdirectory

To define the same commands also for mathjax, follow http://stackoverflow.com/questions/40270302 and add something to `doc/newcommands.js`

adding content
--------------

### chapter

Create file `chapter3.md`, add in the proper location in `$doxyinput` in `Makefile.am`

     input = front.md basics.md bayes.md chapter3.md

Copy over the labeling of sections from `bayes.md`, it doesn't match the doxygen docs on markdown. Add `[TOC]` for a table of contents!

```markdown
Markov chain Monte Carlo {#cha-MCMC}
============

[TOC]

```

### section, subsection

The TOC can get the wrong order if doxygen and markdown syntax for sections and subsections is mixed in one file, so check the TOC carefully. If markdown syntax works for sections, stick with that at all levels!

Section w/o a label don't appear in the TOC as no link can be established, so every section should have a label!

    # A section {#sec-mcmc-a-section}
    ## Factorized proposal {#sec-mcmc-factorized}

### link to a member function

Linking to classes should work directly as in `BCModel` but member functions have to be prefixed by the class as in `BCModel::LogAPrioriProbability`.

### link to section

    @ref sec-mcmc-motiv

### link to next/previous chapter

    @htmlonly
    <br/><br/>
    Return to the <a href="cha-predefined-models.html">previous</a> section or go to the
    <a href="cha-output.html">next</a> section.
    @endhtmlonly

Better would be to define something in javascript that takes two
arguments, the previous and next chapter, and style it with css.

Can't link to installation instructions because its pure markdown for good display on github.

### link to equation

In general not solved http://stackoverflow.com/questions/43191252/references-to-latex-equations-using-doxygen

    \f{align}{
        \label{eq:mc-expect}
        E_P[ f ] = \int \rmdx{ \vecth} P(\vecth) f(\vecth) < \infty .
    \f}

    @latexonly Eq.~\ref{eq:mc-expect-discrete}@endlatexonly

### link to a URL

Contrary to what the doxygen docs at https://www.stack.nl/~dimitri/doxygen/manual/markdown.html#md_links say, a standard markdown link like `[The link text](https://example.net/)` doesn't work. On invoking `doxygen`, there is a warning on `stderr`

    warning: unable to resolve reference to `https:' for \ref command

Just use HTML links

    <a href="https://github.com/bat/bat/blob/master/INSTALL.md">

they are converted properly for latex, too.

### images

Great overview at http://alesnosek.com/blog/2015/06/28/diagrams-and-images-in-doxygen/

From inkscape, export to png for html output with

    inkscape --export-area-drawing --export-png=random-walk.png --export-width 500 random-walk.svg

The images are not created automatically so users don't need to have inkscape etc. installed. Create them with `make histogram.png; make histogram.pdf`.

The markdown image support is not good enough for us, so import images separately for html and latex. Forcibly add both images to git! Also add them to `EXTRA_DIST` in `doc/manual/Makefile.am`

Keep the caption in `""` short because it doesn't look good in HTML and avoid math as mathjax fails in caption. Unfortunately caption has to be copied

    @anchor random-walk-2D
    @image html random-walk.png "2D random walk."
    @image latex random-walk.pdf "2D random walk." width=0.5\textwidth

Refer to the image with

    @ref random-walk-2D "2D example plot"

To set the width of an image in html, use a custom command with three arguments. Omit the third for no caption

    @imageSize{bat.svg,width:300px;,}
    @imageSize{bat.svg,width:300px;,caption}

### copying from introduction.tex

use pandoc to convert to markdown

    pandoc introduction.tex -t markdown -o introduction.md

It's left for emacs regex to translate math expressions.
Display math

     \$\$\\begin{\(.+?\)}\([[:ascii:]]+?\)\\end{.+?}\$\$ -> \\f[
    \1
    \\f]
    \$\$\\begin{\(.+?\)}\([[:ascii:]]+?\)\\end{.+?}\$\$ -> ^J\\f{\1}{\2^J\\f}^J) # display math

Inline math

        \$\([[:ascii:]]+?\)\$ -> \\f$\1\\f$

pandoc uses aligned everywhere and doesn't respect the math environment chosen, so check if outlooks ok, else reinstate what was there before

    {aligned} -> {eqnarray}

Unlabelled section can stay as they are. Convert subsections that have a label.

    ^### \(.+?\) {#\(.+?\)}^J-*$ -> @section \2 \1

Labels may not contain `:` for doxygen

    ^@section section:\(.+? \)\(.+?\)$ -> @section sec-mtf-\1 \2


output specific code
----------

    @htmlonly
    some text only visible in the html output
    @endhtmlonly

    @latexonly
    some text only visible in the latex output
    @endlatexonly

Latex output is copied directly, so `$..$` may not be necessary

### source code examples

Syntax highlighting (and linking in html) with

    @code{.cpp}
    double GetLogPrior(double x)
    {
        return -0.5 * (x - fMean) * (x - fMean) / fSigma / fSigma - log(fSigma) - 0.5 * log(2 * M_PI);
    }
    @endcode

debugging doxygen
-----------------

    doxygen Doxyfile -d

To see warnings better, split output to stdout and stderr

    make > build.log

changing the output style
------------------------

https://www.stack.nl/~dimitri/doxygen/manual/customize.html

### latex

Modify `bat.sty` which is loaded in the preamble.

Todo: `### foo` and `#### bar` look the same, they shouldn't. latex `paragraph` -> `####`

### html

Change names of tabs in `DoxygenLayout.xml`, CSS in
`customdoxygen.css`. We use the twitter bootstrap layout from
https://github.com/Velron/doxygen-bootstrapped defined in javascript
in `doxy-boot.js`.

Warning: a few markups like `@see` and the bibliography are only
displayed in pdf not in html. I still need to figure out how to make
this work with bootstrap. Seems to be a general problem
https://github.com/Velron/doxygen-bootstrapped/issues/21
