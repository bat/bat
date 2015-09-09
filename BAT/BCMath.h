#ifndef __BCMATH__H
#define __BCMATH__H

/*!
 * \namespace BCMath
 * \brief Some useful mathematic functions.
 * \author Frederik Beaujean
 * \author Daniel Greenwald
 * \author Daniel Kollar
 * \author Kevin Kr&ouml;ninger
 * \author Jing Liu
 * \version 1.0
 * \date 08.2008
 * \detail A namespace which encapsulates some mathematical functions
 * necessary for BAT.
 */

/*
 * Copyright (C) 2007-2015, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

// ---------------------------------------------------------
#include <stdexcept>
#include <vector>

class TRandom;

namespace BCMath
{

/** \name Functions for log likelihoods **/
/** @{ */

/**
 * Calculate the natural logarithm of a normal distribution function.
 * @param x point to be evaluated.
 * @param mean Mean.
 * @param sigma Standard deviation.
 * @param norm flag for including normalization constant 1/sqrt(2*pi)/sigma */
double LogGaus(double x, double mean = 0, double sigma = 1, bool norm = false);

/**
 * Calculate the natural logarithm of a normal distribution function with different variances below and above the mode.
 * @param x point to be evaluated.
 * @param mode Mode of the function.
 * @param sigma_below Standard deviation below mode.
 * @param sigma_above Standard deviation above mode.
 * @param norm flag for including normalization constant sqrt(2/pi)/(sigma_below+sigma_above) */
double LogSplitGaus(double x, double mode, double sigma_below, double sigma_above, bool norm = false);

/**
 * Calculate the natural logarithm of a poisson distribution.
 * @param x number of occurances.
 * @param lambda expected number of occurances. */
double LogPoisson(double x, double lambda);

/**
 * Calculates natural logarithm of the Binomial probability using
 * approximations for factorial calculations if calculation for
 * number greater than 20 required using the BCMath::ApproxLogFact
 * function.
 * @param n number of trials.
 * @param k number of successes.
 * @param p probability of success in a trial. */
double LogApproxBinomial(unsigned n, unsigned k, double p);

/**
 * Calculates the logarithm of the nonrelativistic Breit-Wigner distribution. */
double LogBreitWignerNonRel(double x, double mean, double Gamma, bool norm = false);

/**
 * Calculates the logarithm of the relativistic Breit-Wigner distribution. */
double LogBreitWignerRel(double x, double mean, double Gamma);

/**
 * Calculates the logarithm of chi square function:
 * chi2(double x; size_t n)
 */
double LogChi2(double x, int n);

/**
 * Calculates the logarithm of normalized voigtian function:
 * voigtian(double x, double sigma, double gamma)
 *
 * voigtian is a convolution of the following two functions:
 *  gaussian(x) = 1/(sqrt(2*pi)*sigma) * exp(x*x/(2*sigma*sigma)
 *    and
 *  lorentz(x) = (1/pi)*(gamma/2) / (x*x + (gamma/2)*(gamma/2))
 *
 * it is singly peaked at x=0.
 * The width of the peak is decided by sigma and gamma,
 * so they should be positive.
 */
double LogVoigtian(double x, double sigma, double gamma);

/**
* Returns the log of the Gamma PDF.
*/
double LogGammaPDF(double x, double alpha, double beta);

/**
* Return the log of the log normal distribution
*/
double LogLogNormal(double x, double mean = 0, double sigma = 1);

/** @} */

/**
 * Calculates Binomial probability using approximations for
 * factorial calculations if calculation for number greater than 20
 * required using the BCMath::ApproxLogFact function.
 * @param n number of trials.
 * @param k number of successes.
 * @param p probability of success in a trial. */
double ApproxBinomial(unsigned n, unsigned k, double p);

/**
 * Calculates natural logarithm of the Binomial factor "n over k"
 * using approximations for factorial calculations if calculation
 * for number greater than 20 required using the
 * BCMath::ApproxLogFact function.  Even for large numbers the
 * calculation is performed precisely, if n-k < 5
 * @param n upper value in binomial coefficient
 * @param k lower value in binomial coefficient */
double LogBinomFactor(unsigned n, unsigned k);

/**
 * Calculates natural logarithm of the n-factorial (n!) using
 * Srinivasa Ramanujan approximation
 * log(n!) = n*log(n) - n + log(n*(1.+4.*n*(1.+2.*n)))/6. + log(PI)/2.
 * if n > 20.  If n <= 20 it uses BCMath::LogFact to calculate it
 * exactly. */
double ApproxLogFact(double x);

/**
 * Calculates natural logarithm of the Binomial factor "n over k". */
double LogBinomFactorExact(unsigned n, unsigned k);

/**
 * Calculates natural logarithm of the n-factorial (n!) */
double LogFact(unsigned n);

/** Cache factorials for first \arg \c n integers. */
unsigned CacheFactorials(unsigned int n);

/**
 * Returns the nearest integer of a double number. */
int Nint(double x);

/** \name p value methods */
/** @{ */

/**
 * Correct a p value by transforming to a chi^2 with dof=nobservations,
 * then back to a pvalue with dof reduced by number of fit parameters.
 * @param pvalue The p value to correct.
 * @param npar The number of fit parameters.
 * @param nobservations The number of data points.
 * @return corrected p value
 */
double CorrectPValue(const double& pvalue, const unsigned& npar, const unsigned& nobservations) throw (std::domain_error);

/**
 * Calculate the p value using fast MCMC for a histogram and the likelihood as test statistic.
 * The method is explained in the appendix of http://arxiv.org/abs/1011.1674
 *
 * @param observed The counts of observed events
 * @param expected The expected number of events (Poisson means)
 * @param nIterations Controls number of pseudo data sets
 * @param seed Set to nonzero value for reproducible results.
 * @return The p value
 */
double FastPValue(const std::vector<unsigned>& observed, const std::vector<double>& expected,
                  unsigned nIterations = 1e5, unsigned seed = 0) throw (std::invalid_argument);

/** @} */

/** \name Random number generation */
/** @{ */
namespace random
{
/**
 * Chi2 random variate.
 * @param rng Random number generator.
 * @param dof Degree of freedom.
 */
double Chi2(TRandom* rng, double dof);
}
/** @} */
}

// ---------------------------------------------------------

#endif
