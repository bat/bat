Markov chain Monte Carlo {#cha-MCMC}
============

[TOC]

@todo workhorse of BAT. mention Metropolis algorithm, role of proposal function, different choices in bat, why it has to be adapted in prerun. What else happens in prerun:R value and convergence checking. Mention that samples are correlated, chains can get stuck if model wrong/poor.

@section sec-mcmc-motivation Motivation

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

The Markov property implies that the proposal may only on the current
point \f$x\f$ and not on previous points. If the value of \f$\xi\f$ is set
based on a past sequence of iterations of the chain, we need two
stages of sampling in BAT, the *prerun* and the \emph{main
  run}. In the prerun, the chain is run and periodically \f$\xi\f$ is
updated based on the past iterations. In contrast, \f$\xi\f$ is kept fixed
in the main run to have a proper Markov chain. @todo A flow diagram
  might help

BAT offers two kinds of proposal function termed *factorized* and *multivariate*.

@subsection sec-factorized Factorized proposal

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

@subsection sec-multivariate Multivariate proposal

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

@subsection sec-mcmc-prerun Prerun

@todo Explain main options to adjust length of prerun, updates etc.,
discuss proposal-specific options below

The prerun terminates after a maximum number of iterations or if
multiple chains mix according to the R value. Optionally the
minimum/maximum number of iterations can be set with
BCEngineMCMC::SetNIterationsPreRunMin and
BCEngineMCMC::SetNIterationsPreRunMax.

@subsection sec-mcmc-proposal-comparison Comparison

Comparing the factorized proposal to the multivariate proposal, we
generally recommend the multivariate for most purposes.

Use the factorized proposal if you can speed up the computation of the
posterior if you know that some parameters did not change. This can be
useful if the computation is expensive if some but not all
parameters change.
