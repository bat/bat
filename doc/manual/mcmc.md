Markov chain Monte Carlo {#cha-MCMC}
============

[TOC]

@todo workhorse of BAT. mention Metropolis algorithm, role of proposal function, different choices in bat, why it has to be adapted in prerun. What else happens in prerun:R value and convergence checking. Mention that samples are correlated, chains can get stuck if model wrong/poor.

@section sec-mcmc-motiv Motivation

The reason that BAT exists is that nearly any Bayesian analysis
these days is too complicated to be handled analytically. To address
typical questions like

* What is known about a single parameter taking into account
  the uncertainty on all other parameters?
* How are parameters correlated?

one needs to be able to compute and visualize 1D and 2D *marginal distributions*; cf. @ref  sec-marginalization. These are defined as integrals over the posterior; for  example in 2D
\f{align}{
  \label{eq:mcmc-marginal}
  P(\theta_1, \theta_2 \cond D) = \int \prod_{i \ne 1,2} \rmdx{\theta_i} P(\vecth \cond D).
\f}
When the number of parameters grows, the only feasible algorithms to
perform the integration are Monte Carlo methods; i.e., methods based
on random numbers.

The key ingredient in BAT is an implementation of the Metropolis
algorithm to create a Markov chain; i.e. a sequence of (correlated)
samples from the posterior. We use the shorthand MCMC for Markov chain
Monte Carlo.

@section sec-mcmc-foundations Foundations

Efficient MCMC algorithms are the topic of past and current
research. This section is a concise overview of the general idea and
the algorithms available in BAT. For a broader overview, we refer the
reader to the abundant literature; e.g., @todo Robert+Casella, MCMC
  handbook.

In BAT, there are several variants of the random-walk Metropolis Hastings algorithm available. The basic idea is captured in the @ref random-walk-2D "2D example plot". Given an initial point \f$x_0\f$, the Metropolis algorithm produces a sample in each iteration \f$t=1 \dots\f$ as follows:

* Propose a new point \f$y\f$
* Generate a number u from the uniform distribution on [0,1]
* Set \f$x_{t} = y\f$ if \f$ u < \frac{P(y \cond D)}{P(x_{t-1} \cond D)}\f$
* Else stay, \f$x_t = x_{t-1}\f$

@anchor random-walk-2D
@image html random-walk.png "2D random walk with the Metropolis algorithm."
@image latex random-walk.pdf "2D random walk with the Metropolis algorithm." width=0.5\textwidth

In the example plot, the chain begins in the lower left
corner. Rejected moves are indicated by the dashed arrow, accepted
moves are indicated by the solid arrow. The circled number is the
number of iterations the chain stays at a given point \f$\vecth =
(\theta_1, \theta_2)\f$.

@section sec-convergence Convergence
@todo mixing, burn in, R value, multimodal problems

@section sec-mcmc-impl Implementation in BAT

Implementing the Metropolis algorithm, one has to decide on how to
propose a new point based on the current point, that is one needs the
*proposal distribution* \f$q(y \cond x, \xi)\f$ with adjustable parameters
\f$\xi\f$. The main difference between MCMC algorithms is typically given
by different choices of \f$q\f$. The Metropolis algorithm doesn't specify
which \f$q\f$ to choose, so we can and have to select a function \f$q\f$ and
tune \f$\xi\f$ according to our needs.

In BAT, the proposal is *symmetric* around the current point
\f{align}{
  q(y \cond x, \xi) = q(x \cond y, \xi).
\f}

The Markov property implies that the proposal may only depend on the
current point \f$x\f$ and not on previous points. If the value of
\f$\xi\f$ is set based on a past sequence of iterations of the chain,
we need two stages of sampling in BAT, the *prerun* and the *main
run*. In the prerun, the chain is run and periodically \f$\xi\f$ is
updated based on the past iterations. In contrast, \f$\xi\f$ is kept
fixed in the main run to have a proper Markov chain.

## Prerun {#sec-mcmc-prerun}

During the prerun, the proposal is updated. BAT considers three criteria to decide when to end the prerun:

### Efficiency {#sec-mcmc-eff}

The *efficiency*, or acceptance rate, is the ratio  of the accepted over the total number proposal moves. A small efficiency means the chain rarely moves but may then make a large move. A large efficiency means the chain explores well locally but may take a long time to explore the entire region of high probability. Optimality results exists only for very special cases: Roberts and Rosenthal showed that for a Gaussian target with \f$d\f$ independent components and a Gaussian proposal, the optimal target efficiency is 23.4 % for \f$d \geq 5\f$ is  but should be larger for small \f$d\f$; e.g., 44 % is best in one dimension @cite rosenthal2011optimal . Based on our experience, we consider a suitable range of the efficiency to be in \f$[0.15, 0.35]\f$.

@see BCEngineMCMC::SetMinimumEfficiency, BCEngineMCMC::SetMaximumEfficiency

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

@see BCEngineMCMC::SetRValueParametersCriterion Set that maximum R value for all parameters. **Default: 1.1**

@see BCEngineMCMC::SetCorrectRValueForSamplingVariability The strict definition of \f$R\f$ corrects the sampling variability due finite batch size. **Default: false**

@see BCEngineMCMC::GetRValueParameters  \f$R\f$ values are computed during the prerun and they can be retrieved but not set.

### Prerun length {#sec-mcmc-prerun-length}

@see BCEngineMCMC::SetNIterationsPreRunCheck sets the number of iterations between checks


If desired, the statistics can be cleared to remove the effect of a bad initial point with BCEngineMCMC::SetPreRunCheckClear.


@todo A flow diagram might help

BAT offers two kinds of proposal function termed *factorized* and *multivariate*.

## Factorized proposal {#sec-factorized}

@since Factorized was the default and only choice prior to v1.0 and continues to be available

@todo how to activate it

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

This means the posterior is called \f$d\f$ times in every iteration. Since
the acceptance rate

@todo to define before

is typically different from
zero or one, the factorized proposal typically generates a new point
in every iteration that differs from the previous point in some but
not all dimensions.

## Multivariate proposal {#sec-multivariate}

@since 1.0

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

## Prerun {#sec-mcmc-prerun}

@todo Explain main options to adjust length of prerun, updates etc.,
discuss proposal-specific options below

The prerun terminates after a maximum number of iterations or if
multiple chains mix according to the R value. Optionally the
minimum/maximum number of iterations can be set with
BCEngineMCMC::SetNIterationsPreRunMin and
BCEngineMCMC::SetNIterationsPreRunMax.

## Comparison {#sec-mcmc-proposal-comparison}

Comparing the factorized proposal to the multivariate proposal, we
generally recommend the multivariate for most purposes.

Use the factorized proposal if you can speed up the computation of the
posterior if you know that some parameters did not change. This can be
useful if the computation is expensive if some but not all
parameters change.
