Integration {#cha-integration}
============

[TOC]

@section sec-int-motiv Motivation

The reason that BAT exists is that nearly any Bayesian analysis
these days is too complicated to be handled analytically. To address
typical questions like

* What is known about a single parameter taking into account
  the uncertainty on all other parameters?
* How are parameters correlated?

one needs to be able to compute and visualize 1D and 2D *marginal distributions*; cf. @ref  sec-marginalization. These are defined as integrals over the posterior; for  example in 2D
\f{align}{
  \label{eq:int-marginal}
  P(\theta_1, \theta_2 \cond D) = \int \prod_{i \ne 1,2} \rmdx{\theta_i} P(\vecth \cond D).
\f}
To do model comparison, one has to compute the evidence, that is the integral over all parameters
\f{equation}{
Z = \int \rmdx{ \vecth} P(D|\vecth, M) P(\vecth
\cond M).
\f}
Therefore Bayesian inference in practice requires good integration techniques.

For low-dimensional problems, deterministic integration methods are
usually the fastest and most robust but for higher dimensions, Monte
Carlo techniques are the most efficient tools known. Depending on the
setting, BAT defaults to deterministic methods for \f$d \le 2,3\f$
dimensions and uses Monte Carlo for \f$d>3\f$.

@section sec-int-marg Marginalization

The main use case for BAT is to estimate and visualize the marginal
distributions of the posterior. Given a model `m`, all marginal
distributions are estimated as histograms by

@code{.cpp}
m.MarginalizeAll();
@endcode

The individual distributions can be accessed using

@todo snippet of 1D, par 2 vs 5
@todo show how to constrain which distributions are stored, useful if many nuisance parameters.

The method for marginalizing can be selected as follows

@code{.cpp}
m.SetMarginalizationMethod(method);
@endcode

where `method` is an `enum` in  `BCIntegrate::BCMarginalizationMethod`

@todo links to ref guide don't work

`method`  | Details
------------- | -------------
`kMargMetropolis`  | Metropolis algorithm; see @ref cha-mcmc.
`kMargMonteCarlo`  | Sample mean integration in each histogram bin. Least efficient method.
`kMargGrid` | Evaluate target at each bin center. Most efficient for 1D and 2D.
`kMargDefault` | Use `kMargGrid` for \f$d \leq 2 \f$, else `kMargMetropolis`.


@section sec-int-evidence Evidence

In BAT terminology, the evidence or normalization constant of a model
`m` with the default method is computed as

@code{.cpp}
double evidence = m.Normalize();
@endcode

Choose the method of integration explicitly with

@code{.cpp}
double evidence = m.Integrate(method);
@endcode

where `method` is an `enum` in `BCIntegrate::BCIntegrationMethod`.


`BCIntegrationMethod`  | Details
------------- | -------------
`kIntMonteCarlo` | Sample mean. Usually least efficient.
`kIntCuba` | Use a method from cuba.
`kIntGrid` | Approximate Riemann sum over a dense grid for \f$d \leq 3 \f$.
`kIntLaplace` | Laplace approxation. Only approximately valid for peaked unimodal integrands. Usually very fast.
`kIntDefault` | Use Cuba if available. Else use `kIntGrid` for \f$d \leq 3 \f$ and `kIntMonteCarlo` for \f$d > 3 \f$.

The Cuba package itself has four different integration methods. Select the Cuba method with
@code{.cpp}
m.SetCubaIntegrationMethod(method);
m.Integrate();
// alternative
m.IntegrateCuba(method);
@endcode
where `method` is an `enum` in `BCIntegrate::BCCubaMethod`

`BCCubaMethod`  | Details
------------- | -------------
`kCubaVegas` | VEGAS algorithm by Lepage
`kCubaSuave` | Suave algorithm.
`kCubaDivonne` | Divonne algorithm.
`kCubaCuhre` | Cuhre algorithm. Essentially quadrature in higher dimensions. Suffers from curse of dimensionality but most efficient and robust in low dimensions.
`kCubaDefault` | For \f$d=1\f$, use VEGAS, for  \f$d=2,3\f$, use Cuhre, and for  \f$d> 3\f$, use Divonne.

Each Cuba method comes with various parameter values to set. We have
taken over default values from the example that comes with Cuba but
they are by no means optimal for every problem. Please experiment and
consult the Cuba manual that comes with source code on
http://www.feynarts.de/cuba/. All options are accessible in the
namespace `BCCubaOptions`. An example with bogus values

@code{.cpp}
BCCubaOptions::Suave o = m.GetCubaSuaveOptions();
o.flatness = 5;
o.nnew = 5000;
o.nmin = 15;
m.SetNIterationsMax(1e7);
m.SetCubaOptions(o);
m.IntegrateCuba(BCIntegrate::kCubaSuave);
@endcode
