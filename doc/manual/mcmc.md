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

Since samples are not independent, the initial point has some effect on Markov
chain output. The asymptotic results guarantee that, under certain conditions
(see @cite Casella:2004 or @cite brooks2011handbook) a chain of infinite length
is independent of the initial point. In practice, we can only generate a finite
number of points so a decision has to be made when the chain has run long
enough. One helpful criterion is to run multiple chains from different initial
positions and to declare convergence if the chains mixed; i.e. explore the same
region of parameter space. Then the chains have forgotten their initial point.

Non-convergence is a problem that can have many causes including simple bugs in
implementing the posterior. But there are properly implemented posteriors for
which a Markov chain has difficulties to explore the parameter space
efficiently, for example because of strong correlation, degeneracies, or
multiple well separated modes.

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

We implement the adaptive algorithm by Haario et al. @cite Haario:2001,
@cite Wraith:2009if. In brief, the proposal is a multivariate Gaussian or Student's t
distribution whose covariance is learned from the covariance of samples in the
prerun. An overall scale factor is tuned to force the acceptance rate into a
certain range.

the multivariate normal distribution
\f{equation}{
  \label{eq:multivar-normal}
  \mathcal{N}(\vecth | \vecmu, \matsig) = \frac{1}{(2\pi)^{d/2}} \left|\matsig\right|^{-1/2}
  \exp(-\frac{1}{2} (\vecth - \vecmu)^T \matsig^{-1} (\vecth - \vecmu) )
\f}

or the multivariate Student's t distribution
\f{equation}{
  \label{eq:multivar-student}
  \mathcal{T}
  (\vecth | \vecmu, \matsig, \nu)  =
  \frac{\Gamma( (\nu + d) / 2 )}{\Gamma(\nu / 2) (\pi \nu)^{d/2}} \left|\matsig\right|^{-1/2}
  ( 1 + \frac{1}{\nu}(\vecth - \vecmu)^T \matsig^{-1} (\vecth - \vecmu) )^{-(\nu + d)/2}
\f}
can adapt in such a way as to efficiently generate samples from essentially any
smooth, unimodal distribution. The parameter \f$\nu\f$, the degree of freedom,
controls the ``fatness'' of the tails of \f$\mathcal{T}\f$; the covariance of
\f$\mathcal{T}\f$ is related to the scale matrix \f$\matsig\f$ as \f$\frac{\nu}{\nu - 2}
\times \matsig\f$ for \f$\nu > 2\f$, while \f$\matsig\f$ is the covariance of
\f$\mathcal{N}\f$. Hence for finite \f$\nu\f$, \f$\mathcal{T}\f$ has fatter tails than
\f$\mathcal{N}\f$, and for \f$\nu \to \infty\f$, \f$\mathcal{T}(\vecth | \vecmu, \matsig,
\nu) \to \mathcal{N}(\vecth | \vecmu, \matsig)\f$.

Before delving into the details, let us clarify at least qualitatively what we
mean by an efficient proposal.  Our requirements are

* that it allow to sample from the entire target support in finite time,
* that it resolve small and large scale features of the target,
* and that it lead to a Markov chain quickly reaching the asymptotic regime.

An important characteristic of Markov chains is the acceptance rate
\f$\alpha\f$, the ratio of accepted proposal points versus the total length of the
chain. We argue that there exists an optimal \f$\alpha\f$ for a given target and
proposal. If \f$\alpha = 0\f$, the chain is stuck and does not explore the
state space at all. On the contrary, suppose \f$\alpha = 1\f$ and the target
distribution is not globally uniform, then the chain explores only a tiny volume
where the target distribution changes very little. So for some \f$\alpha \in (0,1)\f$,
the chains explore the state space well.

How should the proposal function be adapted? After a chunk of \f$N_{\rm
update}\f$ iterations, we change two things. First, in order to propose points
according to the correlation present in the target density, the proposal scale
matrix \f$\matsig\f$ is updated based on the sample covariance of the last
\f$n\f$ iterations. Second, \f$\matsig\f$ is multiplied with a scale factor
\f$c\f$ that governs the range of the proposal. \f$c\f$ is tuned to force the
acceptance rate to lie in a region of \f$0.15 \le \alpha \le 0.35\f$. The
\f$\alpha\f$ range is based on empirical evidence and the following fact: for a
multivariate normal proposal function, the optimal \f$\alpha\f$ for a normal
target density is \f$0.234\f$, and the optimal scale factor is \f$c =
2.38^2/d\f$ as the dimensionality \f$d\f$ approaches \f$\infty\f$ and the chain
is in the stationary regime @cite Roberts:1997 . We fix the proposal after a
certain number of adaptations, and then collect samples for the final inference
step. However, if the Gaussian proposal function is adapted indefinitely, the
Markov property is lost, but the chain and the empirical averages of the
integrals @latexonly represented by Eq.~\ref{eq:mc-expect}@endlatexonly still
converge under mild conditions @cite Haario:2001.

The efficiency can be enhanced significantly with good initial guesses for \f$c\f$
and \f$\matsig\f$. We use a subscript \f$t\f$ to denote the status after \f$t\f$ updates.
It is often possible to extract an estimate of the target covariance by running
a mode finder like MINUIT that yields the covariance matrix at
the mode as a by product of optimization. In the case of a degenerate target
density, MINUIT necessarily fails, as the gradient is not defined. In such
cases, one can still provide an estimate as
\f{equation}{
  \label{eq:sigma-initial}
  \matsig^0 = \diag ( \sigma_1^2, \sigma_2^2, \dots, \sigma_d^2)
\f}
where \f$\sigma_i^2\f$ is the prior variance of the \f$i\f$-th parameter.  The updated value of
\f$\matsig\f$ in step \f$t\f$ is
\f{equation}{
  \label{eq:sigma-update}
  \matsig^t = (1 - a^t) \matsig^{t-1} + a^t \boldsymbol{S}^t
\f}
where \f$\boldsymbol{S}^t\f$ is the sample covariance of the points in chunk \f$t\f$ and
its element in row \f$m\f$ and column \f$n\f$ is computed as
\f{equation}{
  \label{eq:sample-cov}
  (\boldsymbol{S}^t)_{mn} = \frac{1}{N_{\rm update}-1} \sum_{i=(t-1) \cdot N_{\rm update}}^{t \cdot N_{\rm update}}
  ( (\vecth^i)_m -  \widehat{E_P[(\vecth)_m]} ) ((\vecth^i)_n - \widehat{E_P[(\vecth)_n]} )
\f}
The weight \f$a^t = 1/t^{\lambda}, \lambda \in [0,1]\f$ is chosen to make for a
smooth transition from the initial guess to the eventual target covariance, the
implied cooling is needed for the ergodicity of the chain if the proposal is not
fixed at some point @cite Haario:2001. One uses a fixed value of \f$\lambda\f$,
and the particular value has an effect on the efficiency, but the effect is
generally not dramatic; in this work, we set \f$\lambda=0.5\f$
@cite Wraith:2009if.

We adjust the scale factor \f$c\f$ as described in the pseudocode shown below.
The introduction of a minimum and maximum scale factor is a safeguard against
bugs in the implementation. The only example we can think of that would result
in large scale factors is that of sampling from a uniform distribution over a
very large volume. All proposed points would be in the volume, and accepted, so
\f$\alpha \equiv 1\f$, irrespective of \f$c\f$. All other cases that we
encountered where \f$c > c_{max}\f$ hinted at errors in the code that performs
the update of the proposal.

@code{.cpp}
// default values
αmin = 0.15; αmax = 0.35;
cmin = 1e-5; cmax = 100;
β = 1.5;

// single update of the covariance scale factor
if (α > αmax && c < cmax) {
  c *= β * c
} else if (α < αmin && c > cmin) {
  c /= β
}
@endcode

@see `BCEngineMCMC::SetMultivariateCovarianceUpdateLambda`, `BCEngineMCMC::SetMultivariateEpsilon`, `BCEngineMCMC::SetMultivariateScaleMultiplier`

### Factorized proposal {#sec-mcmc-factorized}

@since Factorized was the default and only choice prior to v1.0 and continues to be available using `BCEngineMCMC::SetProposeMultivariate(false)`

The factorized proposal in \f$d\f$ dimensions is a product of 1D proposals.

We sequentially vary one parameter at a time and complete
one iteration of the chain once a new point has been proposed in *every*
direction. This means the chain attempts to perform a sequence of
axis-aligned moves in one iteration.

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

<!-- @todo A flow diagram might help -->

For the user's convenience, multiple settings related to precision of the Markov
chain can be set at once using `BCEngineMCMC::SetPrecision`. The default setting
is `m.SetPrecision(BCEngineMCMC::kMedium)`.

## Main run {#sec-mcmc-main-run}

In the main run, the proposal is held fixed and each chain is run for `BCEngineMCMC::GetNIterationsRun()` iterations.
@see `BCEngineMCMC::SetNIterationsRun`

To reduce the correlation between samples, a lag can be introduced to take only every 10th element with `BCEngineMCMC::SetNLag`.
