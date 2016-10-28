Sharing / Loading samples {#cha-share}
=================

[TOC]

You can output the samples of BAT's Markov chains in the form of a
ROOT `TTree` saved in a `TFile`. You can read these samples
back into BAT using `BCEmptyModel`:
@code{.cpp}
BCEmptyModel m("MyModel_mcmc.root");
m.Remarginalize();
@endcode
where `MyModel\_mcmc.root` is a `.root` file produced by a BAT
analysis---that is, it contains two `TTree`'s, `X\_mcmc.root`
and `X\_pars.root`, where X is a prefix common to both trees.
BAT
will search the file for two such trees matching the structures
expected from BAT output.
Alternatively, the constructor can be
called with the prefix specified explicitly:
@code{.cpp}
  BCEmptyModel::BCEmptyModel(const std::string& filename, const std::string& name, bool loadObservables).
@endcode
The last argument tells switches on or off the loading of
`BCObservable`'s stored in the `TTree`.

You can also use an instance of your own model instead of `BCEmptyModel` if you provide a
constructor that calls the
@code{.cpp}
BCModel::BCModel(const std::string\&
filename, const std::string\& name, bool loadObservables)
@endcode
constructor.

@todo The links to `BCModel` methods don't work here and elsewhere

Note that instead of calling `MarginalizeAll()`, a reloaded BAT
analysis is marginalized by calling `Remarginalize()`.

All subsequent drawing or accessing of `BCH1D` and `BCH2D`
objects is the same as if the analysis had been run directly rather
than reloaded. For most summary reports of BAT, the actual model
itself is not needed---only the samples. So you can share your results
by sharing your ROOT output files without sharing your model
implementation.

Likewise, when outputting to a ROOT file, BAT (by default) autosaves
the Markov chains to the ouput file at regular intervals. You can use
the steps outlined above to check the state of an ongoing analysis
using `BCEmptyModel`.
