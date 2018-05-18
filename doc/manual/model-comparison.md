Model Comparison {#cha-model-comparison}
================

<!-- @todo I did not mention the `BCDataSet` feature. Should we add it? I think nobody uses it... it was mentioned in the old introduction -->

[TOC]

In case more than one model is defined to explain the same data set the
`BCModelManager` can come in handy. Separate `BCModel`s `m1` and `m2` and their
prior probabilities (70 \% vs 30 \%) are added to the model manager via

@code{.cpp}
BCModelManager mgr;
mgr.AddModel(&m1, 0.7);
mgr.AddModel(&m2, 0.3);
@endcode

For convenience, some of the most important methods for handling `BCModel`s are forwarded to each model in the manager. For example

@code{.cpp}
mgr.SetPrecision(BCEngineMCMC::kQuick);
mgr.MarginalizeAll();
@endcode

is equivalent to

@code{.cpp}
for (unsigned i = 0; i < mgr.GetNModels(); ++i) {
  mgr.GetModel(i)->SetPrecision(BCEngineMCMC::kQuick);
  mgr.GetModel(i)->MarginalizeAll();
}
@endcode

To do model comparison, the evidence is needed for each model, then the Bayes factors and the posterior odds are immediately available. Models are index by an `unsigned`; the first model has index `0` etc.:
@code{.cpp}
mgr.Integrate();
double B = mgr.BayesFactor(0, 1);
mgr.PrintModelComparisonSummary();
@endcode

@see [`examples/advanced/polynomialFit`](https://github.com/bat/bat/tree/master/examples/advanced/advancedGraphFitter)

<!-- @todo If we are fancy, we could include the text output from running the example here. -->
