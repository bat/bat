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

one needs to be able to compute and visualize 1D and 2D *marginal distributions*; cf. @ref  sec-marginalization. These are defined as integrals over the posterior; for  example the 2D marginal distribution is
\f{align}{
  \label{eq:int-marginal}
  P(\theta_1, \theta_2 \cond D) = \int \prod_{i \ne 1,2} \rmdx{\theta_i} P(\vecth \cond D).
\f}
To do model comparison, one has to compute the evidence, that is the integral over all parameters
\f{equation}{
Z = \int \rmdx{ \vecth} P(D|\vecth, M) P(\vecth
\cond M).
\f}
Therefore Bayesian inference in practice requires good integration techniques. These methods are bundled in the class `BCIntegrate` with the exception of the Markov chain code in `BCEngineMCMC`; see also the @ref cha-code-structure.

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

The individual distributions can be accessed using `unsigned` or `std::string` as keys

@code{.cpp}
// one-dimensional
BCH1 d0 = m.GetMarginalized(0);
BCH1D d1 = m.GetMarginalized("name");

// two-dimensional
BCH2D d2 = m.GetMarginalized(3, 4);
BCH2D d3 = m.GetMarginalized("name_of_first_parameter", "name_of_second_parameter");
@endcode

Directly access the underlying histogram with

@code{.cpp}
// one-dimensional
BCH1 h0 = m.GetMarginalizedHistogram(0);

// two-dimensional
BCH2D h1 = m.GetMarginalizedHistogram(3, 4);
@endcode

In case of many (nuisance) parameters, it may be useful to constrain which
histograms are stored because the number of 2D marginals grows like \f$n \cdot
(n-1)/2 \f$ with the number of parameter \f$n\f$. This can be achieved either for all 1D or 2D marginal distributions or for individual combinations. For example:

@code{.cpp}
m.AddParameter("x", 0, 1);
m.AddParameter("y", 0, 1);
m.SetFlagFillHistograms(false); // no 1D marginals for `x`, `y` stored
m.AddParameter("z", 0, 1);
m.SetFillHistogramParPar(0, 2, false); // no `x` vs `z` 2D marginal
@endcode

@see `BCEngineMCMC::SetFillHistogramParPar`, `BCEngineMCMC::SetFillHistogramParObs`, `BCEngineMCMC::SetFillHistogramObsObs`, `BCEngineMCMC::SetFlagFillHistograms`

The method for marginalizing can be selected as follows

@code{.cpp}
m.SetMarginalizationMethod(method);
@endcode

where `method` is an `enum` in  `BCIntegrate::BCMarginalizationMethod`

`method`  | Details
------------- | -------------
`kMargMetropolis`  | Metropolis algorithm; see @ref cha-mcmc.
`kMargMonteCarlo`  | Sample mean integration in each histogram bin. Least efficient method.
`kMargGrid` | Evaluate target at each bin center. Most efficient for 1D and 2D.
`kMargDefault` | Use `kMargGrid` for \f$d \leq 2 \f$, else `kMargMetropolis`.

The availability of marginalization methods can be queried at runtime using `BCIntegrate::CheckMarginalizationAvailability`.

In case any pre- or postprocessing needs to happen to set up data structures, we provide the hooks `virtual void BCIntegrate::MarginalizePreprocess` and `virtual void BCIntegrate::MarginalizePostprocess`. They are empty by default and can be overloaded in a user model.

@section sec-int-evidence Evidence

In BAT terminology, the evidence or normalization constant of a model
`m` with the default method is computed as

@code{.cpp}
double evidence = m.Normalize();
@endcode

Once the normalization has been determined, `BCModel::LogProbability` returns the normalized value of the posterior as opposed to the not normalized result from `BCModel::LogProbabilityNN`.

@note Internally only `BCModel::LogProbabilityNN` is called when integrating or optimizing because it only requires the user to override `BCModel::LogLikelihood` and to set a prior.

`BCModel::Normalize` returns the evidence on the linear scale. Choose the method of integration
explicitly with

@code{.cpp}
double evidence = m.Integrate(method);
@endcode

where `method` is an `enum` in `BCIntegrate::BCIntegrationMethod`.


`BCIntegrationMethod`  | Details
------------- | -------------
`kIntMonteCarlo` | Sample mean. Usually least efficient.
`kIntCuba` | Use a method from cuba.
`kIntGrid` | Approximate Riemann sum over a dense grid for \f$d \leq 3 \f$.
`kIntLaplace` | Laplace approxation. Only approximately valid for peaked unimodal integrands. Incorrect if distribution is heavy tailed or if the mode is near a boundary. Very fast. No uncertainty estimate. Only method implemented on the log scale.
`kIntDefault` | Use Cuba if available. Else use `kIntGrid` for \f$d \leq 3 \f$ and `kIntMonteCarlo` for \f$d > 3 \f$.

General termination criteria for all integration methods except `kIntLaplace` are the desired absolute precision \f$\epsilon_a\f$ and relative precision \f$\epsilon_r\f$ and the minimum and maximum number of iterations:

@code{.cpp}
m.SetAbsolutePrecision(1e-6);
m.SetRelativePrecision(1e-8);
m.SetNIterationsMin(2000);
m.SetNIterationsMax(50000);
@endcode

The integration terminates if
\f{equation}{
|Z-\hat{Z}| \leq \max(\epsilon_a, \epsilon_r Z)
\f}
where \f$\hat{Z}\f$ is the current estimate of the evidence.

@subsection sec-int-rescale Rescaling

In case the log likelihood is very small or very large, going to the linear
scale may take it to exactly zero or infinity in finite precision. This often
happens in practice if the log likelihood is a sum of `N` terms. For example,
assume each factor is 0.5. Then \f$\exp(0.5 N)=\infty\f$ on the computer for
\f$N \geq 1420\f$ using double precision. To avoid this problem, you can rescale
the log likelihood by manually subtracting the value at the mode inside `LogLikelihood`.

Another option, if requirements are satisfied, is to use the Laplace method. It
is naturally implemented on the log scale. While the standard interface to all
integration methods via `BCIntegrate::Normalize()` always transforms to the
linear scale, calling `BCIntegrate::IntegrateLaplace()` directly returns the
evidence on the log scale.

@subsection sec-int-cuba Cuba

The Cuba package itself has four different integration methods. Cuba is an external dependency; see the installation instructions on how to build BAT with Cuba support. If Cuba is available, select the Cuba method with
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

@note Cuba evaluates the posterior in parallel by default. If the posterior is
not thread-safe (see @ref cha-threading), it is recommended to set the
environment variable `CUBACORES` to 1.

Each Cuba method comes with various parameter values to set. We have
taken over default values from the example that comes with Cuba but
they are by no means optimal for every problem. Please experiment and
consult the Cuba manual that comes with the Cuba source code from
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

@section sec-int-slice Slices

To get a quick visualization of a complicated posterior, a slice, or projection, may be preferable to a full marginalization. That is, all parameters except one  [two] are held fixed (instead of integrated over as in marginalization) and the remaining parameter[s] are evaluated on a regular grid. In other words, the conditional distribution

\f{equation}{
P(\theta_1 \cond \vecth_{\backslash 1}, D) = \frac{P(\vecth \cond D)} {P(\vecth_{\backslash 1} \cond D)} .
\f}

The slice is returned as one [two] dimensional histogram. This is achieved with the various variants of `BCIntegrate::GetSlice`.

In a posterior where all parameters are independent, the marginal and conditional distributions coincide because

\f{equation}{
P(\theta_1 \cond \vecth_{\backslash 1}, D) = \frac{P(\vecth_{\backslash 1} \cond D) P(\theta_1 \cond D)} {P(\vecth_{\backslash 1} \cond D)}= P(\theta_1 \cond D) .
\f}

@note Independence rarely holds in practice, so marginalization is still useful.

For example, in a Gaussian posterior with independent parameters, we first find the mode and use these parameter values to find the 1D conditional distribution of the first parameter. In this example, the result is just a Gaussian. Here is how to find the result in general, for non-Gaussian or non-independent posteriors:
@code{.cpp}
#include <TH1.h>
#include <TCanvas.h>
...
int main()
{
...
    m.FindMode(m.GetBestFitParameters());

    TCanvas c;
    unsigned nIterations;
    double log_max;
    int nbins = 100;
    bool normalize = true;
    TH1* slice = m.GetSlice(0, nIterations, log_max, m.GetBestFitParameters(), nbins, normalize);
    slice->Draw("HISTSAME");
    c.Print("slice.pdf");
    delete slice;
...
}
@endcode
