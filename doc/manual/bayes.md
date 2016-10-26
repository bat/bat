Bayesian Statistics {#cha-Bayes}
===================
<!-- @page sec:Bayes Bayesian Statistics -->
[TOC]

In this chapter, we give a concise and necessarily incomplete overview
of Bayesian statistics. Our intention is to cover the bare minimum and
introduce those quantities that appear in BAT. As a toolkit, BAT
allows the user complete freedom in modeling. There is vast literature
on the theory, applications, and modeling; for example, consider
@cite jaynes_probability_2003
@cite hartigan_bayes_1983
@cite sivia_data_2006
@cite mackay_information_2003
@cite dagostini_bayesian_2003
@cite kendall_kendalls_2004
@cite gelman2014bayesian
for a comprehensive overview.

@section sec-basic-terminology Basic terminology
<!-- Basic terminology {#sec-basic-terminology} -->
<!-- ## Basic terminology {#sec-basic-terminology} -->

Bayesian statistics is a framework for quantitative inductive or
plausible reasoning; i.e., the optimal processing of incomplete
information. The basic idea is to associate a degree of belief to
logical propositions; e.g., the mass of the Higgs boson is 125
GeV. From fundamental axioms about reasoning, one can then show that
the calculus of degree of belief is simply the ordinary calculus of
probability theory; see @cite jaynes_probability_2003 for a thorough
discussion. Deductive reasoning is included as the limiting case in
which the degree of belief is either 0 or 1.

From now on, we will use the symbol \f$P(A)\f$ to denote both the *degree
of belief* in proposition \f$A\f$ and the *probability* of \f$A\f$. In our
applications below, \f$A\f$ is often one value out of a continuum, so we
use \f$P(A)\f$ also to denote the *probability density* of \f$A\f$.  The
*conditional probability* of \f$A\f$ given \f$B\f$ is \f$P(A \cond B)\f$.

The two central tasks of the natural sciences are to learn about
nature from data and to make predictions for (future) experiments. A
*model* \f$M\f$ is a proxy for all discrete pieces of information
relevant to calculating the degree of belief. The model can contain
*parameters* \f$\vecth\f$ that take on values in a continuum, perhaps
subject to constraints as for example \f$\scath_1 \geq 0\f$. Bayesian
reasoning provides an update rule to adjust the degree of belief based
on new information available in the form of *observed data*
\f$D\f$. This update rule is the celebrated *Bayes' theorem*
\f{align}{ \label{eq:bayes-thm} \boxed{ P(\vecth \cond D, M) \propto
P(D|\vecth, M) P(\vecth \cond M)} \f} \f$P(\vecth \cond M)\f$ is the
*prior density*, \f$P(D|\vecth, M)\f$ is called the *probability of
the data* when treated as a function of \f$D\f$, and known as the
*likelihood* when considering the dependence on \f$\vecth\f$, for
fixed \f$D\f$. The model-dependent normalization constant is known as
the *evidence* or *marginal likelihood*: \f{equation}{
\label{eq:evidence} Z = \int \rmdx{ \vecth} P(D|\vecth, M) P(\vecth
\cond M).  \f} Finally, the left-hand side \f$P(\vecth | D, M)\f$ is
the *posterior density*. Prior and posterior (*density* is usually
omitted) represent the state of knowledge about the parameter
\f$\vecth\f$ before and after seeing the data. Note that \f$\vecth\f$
appears on opposite sides of `|` in \f$P(D|\vecth,M)\f$ and
\f$P(\vecth|D,M)\f$. That's why Bayes' theorem is also known as the
theorem of *inverse probability*.

@subsection sec-marginalization Marginalization

Suppose there are two parameters, \f$\vecth = (\scath_1, \scath_2)\f$, and
\f$\scath_1\f$ is the *parameter of interest* whereas \f$\scath_2\f$ is a
*nuisance parameter*. In Bayes' theorem, there is no fundamental
distinction between parameters of interest and nuisance parameter,
they are all just parameters. But often the goal of the analysis is to
extract the posterior of \f$\scath_1\f$ while \f$\scath_2\f$ is only needed at
an intermediate stage; for example in order to correctly model the
measurement process of \f$D\f$. From the joint posterior \f$P(\scath_1,
\scath_2 | D)\f$, we compute the *marginalized* posterior and can
remove the dependence on \f$\scath_2\f$ by integration
\f{align}{
  \label{eq:marginal}
  P(\scath_1 | D) = \int \rmdx{\scath_2} P(\scath_1, \scath_2 | D).
\f}

@subsection sec-model-comparison Model comparison

If there is only a single model under consideration, and no potential for
cconfusion, the model label \f$M\f$ is implied and usually omitted from the
equations. But suppose that there are two competing models, \f$M_1, M_2\f$, with
parameters \f$\scath_{1,2}\f$, that quantitatively predict the outcome \f$D\f$ of an
experiment. The task is to find the model with the higher degree of
belief. Using Bayes' theorem, the *posterior odds* of the models are easily
found as
\f{equation}{
  \label{eq:post-odds}
  \frac{P(M_1|D)}{P(M_2|D)}
   = B_{12}  \cdot  \frac{P(M_1)}{P(M_2)},
\f}
 where the *Bayes factor* of \f$M_1\f$ versus \f$M_2\f$, \f$B_{12}\f$, is just the ratio
of the evidences
\f{equation}{
  \label{eq:Bayes-factor}
  B_{12}= \dfrac{P(D|M_1)}{P(D|M_2)} = \frac{Z_1}{Z_2}
  = \frac{\int \rmdx{\vecth_1} P(D|\vecth_1, M_1) P(\vecth_1, M_1)}
  {\int \rmdx{ \vecth_2} P(D|\vecth_2, M_2) P(\vecth_2, M_2)}
\f}
The *prior odds* \f$P(M_1)/P(M_2)\f$ represent the relative degree of belief
in the models, independent of the data.
% The data are accounted for in the Bayes factor.
The Bayes factor quantifies the relative shift of degree of belief
induced by the data. In general, \f$\dim \vecth_1 \ne \dim \vecth_2\f$,
and without loss of generality let \f$\dim \vecth_1 < \dim
\vecth_2\f$. The Bayes factor automatically penalizes \f$M_2\f$ for its
larger complexity, as the prior mass is spread out over a
higher-dimensional volume. However, this can be compensated if the
likelihood \f$P(D|\vecth_2, M_2)\f$ is significantly higher in regions of
reasonably high prior density; i.e. the Bayes factor implements
Occam's razor the simplest model that describes the observations is
preferred.

@section sec-goodness-fit Goodness of fit

In the Bayesian
approach, there is, however, no straightforward answer to the
following question: if there is only one model at hand, how to decide
if that model is sufficient to explain the data, or if the search for
a better model needs to continue?  The standard procedure to tackle
this problem of evaluating the *goodness of fit} is to define a
test statistic \f$T=T(D)\f$ and to evaluate the following tail-area
probability, that is the \f$p\f$ value
\f{align}{
  \label{eq:pvalue}
  p \equiv \int_{T>T_{obs}} \rmdx{T} P(T|M).
\f}
Care has to be taken in the usage, computation, and interpretation of
\f$p\f$ values. An introduction a Bayesian interpretation of \f$p\f$ values
with applications in BAT is available at @cite Beaujean:2011zz ; see
also the references therein.

@section sec-representation-bat Representation in BAT

**Fred: the whole section might go somewhere else where users find it more
  easily. I fear it would get lost here**

In BAT, a model \f$M\f$ is represented as a C++ subclass of @ref
BCModel. The crucial parts are to define the likelihood \f$P(D \cond
\vecth, M)\f$ @ref BCModel::LogLikelihood, the prior \f$P(\vecth \cond
M)\f$ BCModel::LogAPrioriProbability, and the parameters
\f$\vecth\f$. From these quantities, BAT can compute the unnormalized
posterior \f$P(\vecth \cond D, M)\f$ BCModel::LogProbabilityNN. To
avoid numerical overflow, BAT operates on the log scale whenever
possible.  The key methods of @ref BCModel that a user has to
implement are

    virtual double LogLikelihood(const std::vector<double>& params)
    virtual double LogAPrioriProbability(const std::vector<double>& params)

<!-- BCEngine link works but BCModel method not found. Why? -->
<!-- @see BCEngineMCMC::SetNIterationsPreRunMax -->
@see BCModel::LogAPrioriProbability BCModel::LogLikelihood BCModel::LogProbabilityNN

The parameter values are passed in simply as numbers to likelihood and
prior, all parameters are assumed to be real and continuous. Discrete
parameters are not supported. The support of \f$\vecth\f$ is a
hyperrectangle whose bounds are given by the bounds of the individual
parameters when added to the model with BCModel::AddParameter

    bool BCModel::AddParameter(const std::string& name, double min, double max,
                               const std::string& latexname = "",
                               const std::string& unitstring = "")

The optional `latexname` and `unitstring` are used only
for labeling plot axes; it's intended usage is to pretty up plots. For
example, a parameter `theta` representing a time measured in
seconds is defined as

    AddParameter("theta", 0, 1, "#theta", "s");

and whenever `theta` appears on the axis of a plot, it will appear as
\f$\theta\f$ [s]. Note that the plots are created with \ROOT, so the
`latexname` has to be in \ROOT syntax which is basically \f$LaTeX\f$
syntax but the backslash `\` is replaced by the hash `#`.
