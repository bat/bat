 <!-- First header is ignored by doxygen for the main page -->


@anchor bat-logo
@imageSize{bat.svg,width:300px;,}
@image latex bat.pdf width=0.3\textwidth

### For the impatient

For a quick tutorial covering the basic features, jump directly to
the [basic tutorial](@ref cha-basics).

The purpose of this document is to introduce the Bayesian Analysis
Toolkit (BAT), a C++ package providing tools to

* compare model predictions with data,
* draw conclusions on the validity of a model as a representation of the data,
* and to extract the values of the free parameters of a model.

BAT was originally developed at
  the [Max Planck institute for physics](http://mpp.mpg.de) in Munich,
  Germany, in the context of particle physics. BAT uses [ROOT] for
  data handling and graphical output.

### How to read this document

If you are new to BAT, we recommend to start with the [basic tutorial](@ref
cha-basics) to quickly get you going. If you then want to learn more, check the
brief introduction to the statistical basics in @ref cha-Bayes. The core
algorithm of BAT is Markov chain Monte Carlo, described in @ref cha-mcmc.
Further chapters describe other tools of BAT, such as @ref cha-optimization,
specific statistical models, and advanced features. <!-- @htmlonly Use the
navigation bar on the left to jump to the chapters @endhtmlonly -->

### Further help

Please visit our [home page][BAThome] for an overview of papers,
research, and other activities related to BAT. The code development
takes place over at [github][BATgithub].

This document is meant to is quickly enable to you to find your way around BAT.
We describe the common workflow and explain how to solve some not-so-common
tasks. To keep the size of this document manageable, we refer to the [html
reference guide](../../ref-guide/html/index.html) in which all classes, methods
etc. are documented. The BAT source code comes with many examples that
illustrate how to use the various features. You can browse them [online](https://github.com/bat/bat/tree/master/examples).

If you run into a problem with BAT, the preferred method is to create a new
issue on [github](https://github.com/bat/bat/issues). This allows other users to
easily find your issue online and benefit from the solution. If you provide some
code for us to reproduce the problem, we can help you much faster. Have a look
at issues with the label `troubleshooting` or `question`.

For general questions or comments that do not belong into the public in the form
of issues, you can contact the developers at bat@mpp.mpg.de.

[BAThome]: http://mpp.mpg.de/bat "BAT homepage"
[BATgithub]: https://github.com/bat/bat "BAT github"
[BATref]: http://mpp.mpg.de/bat/docs/refman/latest/ "BAT reference guide"
[ROOT]: https://root.cern.ch/ "ROOT homepage"
