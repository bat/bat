Defining a factorized prior {#cha-factorized}
================

[TOC]

If a model does not overload `BCModel::LogAPrioriProbability`,
then its prior is the product of the individual priors for each of its
parameters. We call these the factorized priors. To set a parameter's
prior call
@code{.cpp}
BCParameter::SetPrior(BCPrior* const prior);
@endcode
You can call factorized priors within your overloaded
`LogAPrioriProbability` by querying the log of the prior of a
particular parameter with
@code{.cpp}
BCParameter::GetLogPrior(double x)
@endcode
The prior you set for a parameter need only inherit from
`BCPrior`. BAT has several built in prior classes, such as
`BCConstantPrior` and `BCGaussianPrior`.  You can implement
new factorized priors by inheriting from the `BCPrior` class. To
illustrate how to do this, we will work through the construction of
the Gaussian prior:
@code{.cpp}
class BCGaussianPrior : public BCPrior
@endcode
`BCPrior` is a pure-virtual class with three methods that must be overloaded:
@code{.cpp}
double GetLogPrior(double x);
BCPrior* clone() const;
bool IsValid() const;
@endcode
The first one contains the meat of our new prior:
@code{.cpp}
double GetLogPrior(double x)
{
  return -0.5 * (x - fMean) * (x - fMean) / fSigma / fSigma - log(fSigma) - 0.5 * log(2 * M_PI);
}
@endcode
This is, naturally, the log of a Guassian prior. The second one should
simply return a copy of the prior:
@code{.cpp}
BCPrior* clone() const
{
  return new BCGaussianPrior(*this);
}
@endcode
The last function is required for checking that everything that is
needed by the prior is properly set. You'll notice that our example
calls on two member variables in its `GetLogPrior`: `fMean`
and `fSigma`. So we must check whether they have been properly set:
@code{.cpp}
bool IsValid() const
{
  return std::isfinite(fMean) and std::isfinite(fSigma) and fSigma > 0;
}
@endcode
This checks that the mean and standard deviation are both finite and
that the standard deviation is positive semi-definite, as it need be.

Naturally to be of use, we also create a constructor that allows us to
set the mean and standard deviation at creation; and getters and
setters that allow us to access them. (See the source code of
`BCGaussianPrior` for their implementation.)

This is all that is required to create a new prior. `BCPrior` has
several methods for getting properties of the prior: the mode; the
integral over a range; the \f$n\f$'th raw, central, and standardized
moments; the mean; the variance; the standard deviation; the skewness;
and the kurtosis. These functions make their own internal calculations
and need not be overloaded. If you wish to speed them up, or provide
exact results, you may overload them. But take care that you overload
them to give back proper results! For example, the mode of the
Gaussian distribution is not simply `fMean`, since this presumes
the range over which we query contains `fMean`. The more proper
implementation is
@code{.cpp}
double GetMode(double xmin, double xmax)
{
  if (fMean < xmin)
    return xmin;
  if (fMean > xmax)
    return xmax;
  return fMean;
}
@endcode
