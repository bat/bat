<!-- First header is ignored by doxygen for the main page -->
<!-- ### Overview -->

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

### Structure of this document

After the [basic tutorial](@ref cha-basics), we give a brief introduction
to @ref cha-Bayes and describe the core algorithm of BAT in @ref
cha-mcmc. Further chapters describe other tools of BAT, specific
statistical models, and advanced features.
<!-- @htmlonly Use the navigation bar on the left to jump to the chapters @endhtmlonly -->

### Further help

Please visit our [home page][BAThome] for an overview of papers,
research, and other activities related to BAT. The code development
takes place [github][BATgithub].

This document is meant to is quickly enable to you to find your way
around BAT. We describe the common workflow and explain how to solve
some not-so-common tasks. To keep the size of this document
manageable, we refer to the [html reference guide](../../ref-guide/html/index.html) in which
all classes, methods etc. are documented.

If you run into a problem with BAT, the preferred method is to create
a new issue on [github][BATgithub]. This allows other users to easily
find your issue online and benefit from the solution. If you provide
some code for us to reproduce the problem, we can help you much
faster.

For general questions or comments that do not belong into the public
in the form of issues, you can contact the developers at
bat@mpp.mpg.de.

[BAThome]: http://mpp.mpg.de/bat "BAT homepage"
[BATgithub]: https://github.com/bat/bat "BAT github"
[BATref]: http://mpp.mpg.de/bat/docs/refman/latest/ "BAT reference guide"
[ROOT]: https://root.cern.ch/ "ROOT homepage"
