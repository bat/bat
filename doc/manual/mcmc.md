Markov chain Monte Carlo {#cha-mcmc}
=======================

[TOC]

# Motivation {#sec-mcmc-motiv}

Among the integration methods introduced in @ref cha-integration, the Monte Carlo method is the most powerful one in high dimensions. The term Monte Carlo is used as a synonym for the use of pseudo-random numbers. Markov chains are a particular class of Monte Carlo algorithms designed to generate correlated samples from an arbitrary distribution. The central workhorse in BAT is an adaptive Markov chain Monte Carlo (MCMC) implementation based on the Metropolis algorithm. It allows users to marginalize a posterior without requiring manual tuning of algorithm parameters. In complicated cases, tweaking the parameters can substantially increase the efficiency, so BAT gives users full access to all tuning parameters.

# Foundations {#sec-mcmc-foundations}

## Monte Carlo integration {#sec-mcmc-integration}

We begin with the *fundamental Monte Carlo* principle. Suppose we have a posterior
probability density \f$P(\vecth|D)\f$, often called the *target*
density, and an arbitrary function \f$f(\vecth)\f$ with finite
expectation value under \f$P\f$
\f{align}{
    \label{eq:mc-expect}
    E_P[ f ] = \int \rmdx{ \vecth} P(\vecth|D) f(\vecth) < \infty .
\f}
Then a
set of draws \f$\{ \vecth_i:i=1 \dots N \}\f$ from the density \f$P\f$, that is \f$\vecth_i \sim P\f$, is
enough to estimate the expectation value. Specifically, the integral
can be replaced by the estimator (distinguished by the symbol
\f$\widehat{\phantom{a}}\f$)
\f{align}{ \label{eq:mc-expect-discrete} \boxed{
\widehat{E_P[f]} \approx \frac{1}{N} \sum_{i=1}^{N} f(\vecth_i), \;
\vecth \sim P  }
\f}
As \f$N \to \infty\f$, the estimate
converges almost surely at a rate \f$\propto 1/\sqrt{N}\f$
by the strong law of large numbers if \f$\int \rmdx{ \vecth}
P(\vecth|D) f^2(\vecth) < \infty\f$ @cite Casella:2004 . This is true for independent samples from the target but also for correlated samples. The only thing that changes is the increased variance of the estimator due to correlation.

@anchor mcmc-histogram
@imageSize{histogram.svg,width:300px;,Histogram approximation to the 1D marginal.}
@image latex histogram.pdf "Histogram approximation to the 1D marginal." width=0.3\textwidth

How does this @latexonly Eq.~\ref{eq:mc-expect-discrete}@endlatexonly
relate to Bayesian inference? Upon applying Bayes’ theorem to
real-life problems, one usually has to marginalize over several
parameters, and this can usually not be done analytically, hence one
has to resort to numerical techniques. In low dimensions, say \f$d \le
2\f$, quadrature and other grid-based methods are fast and accurate,
but as \f$d\f$ increases, these methods generically suffer from the
*curse of dimensionality*. The number of function evaluations grows
exponentially as \f$\order{m^d}\f$, where \f$m\f$ is the number of
grid points in one dimension. Though less accurate in few dimensions,
Monte Carlo — i.e., random-number based — methods are the first choice
in \f$d \gtrsim 3\f$ because the computational complexity is (at least
in principle) independent of \f$d\f$. Which function \f$f\f$ is of
interest to us? For example when integrating over all but the first
dimension of \f$\vecth\f$, the marginal posterior probability
that \f$\theta_1\f$ is in \f$[a,b)\f$ can be estimated as

\f[
\label{eq:disc-marg} P(a \le \theta_1 \le b|D) \approx \frac{1}{N}
\sum_{i=1}^{N} \mathbf{1}_{\theta_1 \in [a,b)} (\vecth_i) \, \f]
with the *indicator function*
\f[ \label{eq:indicator-fct}
\mathbf{1}_{\theta_1 \in [a,b)} (\vecth) = \begin{cases} 1, \theta_1
\in [a,b) \\ 0, {\rm else} \end{cases}
\f]

This follows immediately from the Monte Carlo principle with
\f$f(\vecth) = \mathbf{1}_{\theta_1 \in [a,b)}(\vecth)\f$. The major
simplification arises as we perform the integral over \f$d-1\f$
dimensions simply by ignoring these dimensions in the indicator
function. If the parameter range of \f$\theta_1\f$ is partitioned into
bins, then the above holds in every bin, and defines the histogram
approximation to \f$P(\theta_1|D)\f$. In exact analogy, the 2D
histogram approximation is computed from the samples for 2D bins in
the indicator function. For understanding and presenting the results
of Bayesian parameter inference, the set of 1D and 2D marginal
distributions is the primary goal. Given samples from the full
posterior, we have immediate access to *all* marginal distributions at
once; i.e., there is no need for separate integration to obtain for
example \f$P(\theta_1|D)\f$ and \f$P(\theta_2|D)\f$. This is a major
benefit of the Monte Carlo method in conducting Bayesian inference.

## Metropolis algorithm {#sec-mcmc-metropolis}

The key ingredient in BAT is an implementation of the Metropolis
algorithm to create a Markov chain; i.e. a sequence of (correlated)
samples from the posterior. We use the shorthand MCMC for Markov chain
Monte Carlo.

Efficient MCMC algorithms are the topic of past and current
research. This section is a concise overview of the general idea and
the algorithms available in BAT. For a broader overview, we refer the
reader to the abundant literature; e.g.,
@cite Casella:2004
@cite brooks2011handbook .

In BAT, there are several variants of the random-walk Metropolis Hastings algorithm available. The basic idea is captured in the @ref random-walk-2D "2D example plot". Given an initial point \f$\vecth_0\f$, the Metropolis algorithm produces a sample in each iteration \f$t=1 \dots N \f$ as follows:

* Propose a new point \f$\tilde{\vecth}\f$
* Generate a number u from the uniform distribution on [0,1]
* Set \f$\vecth_{t} = \tilde{\vecth}\f$ if \f$ u < \frac{P(\tilde{\vecth} \cond D)}{P(\vecth_{t-1} \cond D)}\f$
* Else stay, \f$\vecth_t = \vecth_{t-1}\f$

@anchor random-walk-2D
@image html random-walk.png "2D random walk with the Metropolis algorithm."
@image latex random-walk.pdf "2D random walk with the Metropolis algorithm." width=0.5\textwidth

In the example plot, the chain begins in the lower left
corner. Rejected moves are indicated by the dashed arrow, accepted
moves are indicated by the solid arrow. The circled number is the
number of iterations the chain stays at a given point \f$\vecth =
(\theta_1, \theta_2)\f$.

In each iteration \f$t\f$, one updates the estimate of the 1D marginal
distribution \f$P(\theta_1 | D)\f$ by adding the first coordinate of
\f$\vecth_t\f$ to a histogram. Repeat this for all other coordinates to update
the other \f$(d-1)\f$ 1D marginals. And redo it for all pairs of coordinates to
estimate the 2D marginals.

As a concrete example, suppose the chain has 5 iterations in 2D:

\f$t\f$ | \f$\vecth_t\f$
:----- | :----------
1       | \f$(1.1, 2.3)\f$
2       | \f$(1.1, 2.3)\f$
3       | \f$(3.8, 1.8)\f$
4       | \f$(2.4, 5.2)\f$
5       | \f$(1.8, 4.2)\f$

Let us choose a histogram to approximate the 1D marginal posterior for
\f$\theta_1\f$ with five bins from \f$[n, n+1)\f$ for \f$n=0 \dots 4\f$ such
that the right edge of the bin is not included. Up to \f$t=5\f$, the histogram
is

 \f$n\f$| weight
:----- | :----------
0       | 0
1       | 3
2       | 1
3       | 1
4       | 0

In the end, we usually normalize the histogram so it estimates a proper
probability density that integrates to 1.

## Convergence {#sec-mcmc-convergence}
@todo mixing, burn in, R value, multimodal problems

Since samples are not independent, the initial point has some effect on Markov
chain output. The asymptotic results guarantee that, under certain conditions
(see @cite Casella:2004 or @cite brooks2011handbook) a chain of infinte length
is independent of the initial point. In practice, we can only generate a finite
number of points so a decision has to be made when the chain has run long
enough. One helpful criterion is to run multiple chains from different initial
positions and to declare convergence if the chains mixed; i.e. explore the same
region of parameter space. Then the chains have forgotten their initial point.

# Implementation in BAT {#sec-mcmc-impl}

Implementing the Metropolis algorithm, one has to decide on how to
propose a new point based on the current point, that is one needs the
*proposal function* \f$q(\tilde{\vecth} \cond \vecth_t, \xi)\f$ with adjustable parameters
\f$\xi\f$. The main difference between MCMC algorithms is typically given
by different choices of \f$q\f$. The Metropolis algorithm doesn't specify
which \f$q\f$ to choose, so we can and have to select a function \f$q\f$ and
tune \f$\xi\f$ according to our needs.

In BAT, the proposal is *symmetric* around the current point
\f{align}{
  q(\tilde{\vecth} \cond \vecth_t, \xi) = q(\vecth_t \cond \tilde{\vecth}, \xi).
\f}

The Markov property implies that the proposal may only depend on the current
point \f$\vecth_t\f$ and not on any previous point. If the value of \f$\xi\f$ is
set based on a past sequence of iterations of the chain, we need two stages of
sampling in BAT, the *prerun* and the *main run*. In the prerun, the chain is
run and periodically \f$\xi\f$ is updated based on the past iterations. In
contrast, \f$\xi\f$ is kept fixed in the main run to have a proper Markov chain.

## Proposal functions {#sec-mcmc-proposal}

BAT offers two kinds of proposal function termed *factorized* and
*multivariate*. The general form is either a Gaussian or Student's t
distribution. In the factorized case, the joint distribution is a product of 1D
distributions. In the multivariate case, a dense covariance matrix is used that
allows correlated proposals. In either case, the default is Student's t
distribution with one degree of freedom (dof); i.e., a Cauchy distribution.
Select the proposal like this:

@code{.cpp}
m.SetProposeMultivariate(true);
m.SetProposalFunctionDof(5); // Student's t with 5 degrees of freedom
m.SetProposalFunctionDof(-1); // Gaussian
@endcode

### Multivariate proposal {#sec-mcmc-multivariate}

@since Introduced and set as the default in v1.0

Changing all \f$d\f$ parameters at once within one iteration is an
all-or-nothing approach. If the proposed move is accepted, all
parameters have changed for the price of a single evaluation of the
posterior. If the move is rejected, the new point is identical to the
old point and the chain does not explore the parameter space.

We implement the adaptive algorithm @cite Haario:2001 by Haario et al. In
brief, the proposal is a multivariate Gaussian or Student's t
distribution whose covariance is learned from the covariance of
samples in the prerun. An overall scale factor is tuned to force the
acceptance rate into a certain range.

@todo Pseudocode
@todo copy from my thesis

@see `SetMultivariateCovarianceUpdateLambda`, `SetMultivariateEpsilon`, `SetMultivariateScaleMultiplier`

### Factorized proposal {#sec-mcmc-factorized}

@since Factorized was the default and only choice prior to v1.0 and continues to be available using `BCEngineMCMC::SetProposeMultivariate(false)`

The factorized proposal in \f$d\f$ dimensions is a product of 1D proposals.

We sequentially vary one parameter at a time and complete
one iteration of the chain once a new point has been proposed in *every*
direction. This means the chain attempts to perform a sequence of
axis-aligned moves in one iteration.

@todo Easiest to understand would be pseudocode

Each 1D proposal is a Cauchy or Breit-Wigner function centered on the
current point. The scale parameter is adapted in the prerun to achieve
an acceptance rate in a given range that can be adjusted by the
user. Note that there is a separate scale parameter in every dimension.

This means the posterior is called \f$d\f$ times in every iteration. Since the
acceptance rate is typically different from zero or one, the factorized proposal
typically generates a new point in every iteration that differs from the
previous point in some but not all dimensions.

### Comparison {#sec-mcmc-proposal-comparison}

Comparing the factorized proposal to the multivariate proposal, we
generally recommend the multivariate for most purposes.

Use the factorized proposal if you can speed up the computation of the
posterior if you know that some parameters did not change. This can be
useful if the computation is expensive if some but not all
parameters change.

## Prerun {#sec-mcmc-prerun}

During the prerun, the proposal is updated. BAT considers three criteria to decide when to end the prerun. The prerun takes some minimum number of iterations and is stopped no matter what if the maximum number of prerun iterations is reached. In between, the prerun terminates if the efficiency and the \f$R\f$ value checks are ok. To perform the prerun manually, do

@code{.cpp}
m.MetropolisPreRun();
@endcode

In most cases this is not needed, because `BCModel::MarginalizeAll` calls
`BCEngineMCMC::Metropolis`, which will take care of the prerun and the main run
and all the data handling associated with it. To force a prerun to be run again, perhaps after it failed and some settings have been adjusted, do:

@code{.cpp}
m.SetFlagPreRun(true);
m.MarginalizeAll();
@endcode

### Efficiency {#sec-mcmc-eff}

The *efficiency*, or acceptance rate, is the ratio  of the accepted over the total number proposal moves. A small efficiency means the chain rarely moves but may then make a large move. A large efficiency means the chain explores well locally but may take a long time to explore the entire region of high probability. Optimality results exists only for very special cases: Roberts and Rosenthal showed that for a Gaussian target with \f$d\f$ independent components and a Gaussian proposal, the optimal target efficiency is 23.4 % for \f$d \geq 5\f$ is  but should be larger for small \f$d\f$; e.g., 44 % is best in one dimension @cite rosenthal2011optimal . Based on our experience, we use a default range for the efficiency as \f$[0.15, 0.35]\f$.

@see `BCEngineMCMC::SetMinimumEfficiency`, `BCEngineMCMC::SetMaximumEfficiency`

### R value {#sec-mcmc-Rvalue}

The *R value* @cite Gelman:1992 by Gelman and Rubin quantifies the
estimated scale reduction of the uncertainty of an expectation value
estimated with the samples if the chain were run infinitely
long. Informally, it compares the mean and variance of the expectation
value for a single chain with the corresponding results of multiple
chains. If the chains mix despite different initial values, then we
assume that they are independent of the initial value, the burn-in is
over, and the samples produce reliable estimates of quantities of
interest. For a single chain, the \f$R\f$ value cannot be computed.

In BAT, we monitor the expectation value of each parameter and declare
convergence if all R values are below a threshold. Note that the R
values are estimated from batches of samples, and they usually
decrease with more iterations but they may also increase, which
usually is a clear indication that the chains do not mix, perhaps due
to multiple modes that trap the chains.

@see `BCEngineMCMC::SetRValueParametersCriterion` Set the maximum allowed \f$R\f$ value for all parameters. By defition, \f$R\f$ cannot go below 1 except for numerical inaccuracy. **Default: 1.1**

@see `BCEngineMCMC::SetCorrectRValueForSamplingVariability` The strict definition of \f$R\f$ corrects the sampling variability due finite batch size. **Default: false**

@see `BCEngineMCMC::GetRValueParameters`  \f$R\f$ values are computed during the prerun and they can be retrieved but not set.

### Prerun length {#sec-mcmc-prerun-length}

Defining convergence automatically based on the efficiency or the \f$R\f$ value is convenient may be too conservative if the user knows a good initial value, a good proposal, etc. For more control the minimum and maximum length of the prerun can be set, too.

@see `BCEngineMCMC::SetNIterationsPreRunMin`, `BCEngineMCMC::SetNIterationsPreRunMax`

@see `BCEngineMCMC::SetNIterationsPreRunCheck` sets the number of iterations between checks

If desired, the statistics can be cleared to remove the effect of a bad initial point with `BCEngineMCMC::SetPreRunCheckClear` after some set of iterations

@todo A flow diagram might help

For the user's convenience, multiple settings related to precision of the Markov
chain can be set at once using `BCEngineMCMC::SetPrecision`. The default setting
is `m.SetPrecision(BCEngineMCMC::kMedium)`.

## Main run {#sec-mcmc-main-run}

@see `SetNIterationsRun`
