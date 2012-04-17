#ifndef __BCMATH__H
#define __BCMATH__H

/*!
 * \namespace BCMath
 * \brief Some useful mathematic functions.
 * \author Daniel Kollar
 * \author Kevin Kr&ouml;ninger
 * \author Jing Liu
 * \author Frederik Beaujean
 * \version 1.0
 * \date 08.2008
 * \detail A namespace which encapsulates some mathematical functions
 * necessary for BAT.
 */

/**
 * Copyright (C) 2008-2012, Daniel Kollar and Kevin Kroeninger.
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

// ---------------------------------------------------------

#define BCMATH_NFACT_ALIMIT 20

// ---------------------------------------------------------
#include <stdexcept>
#include <vector>

class TH1D;

namespace BCMath
{

   /**
    * Calculate the natural logarithm of a gaussian function with mean
    * and sigma.  If norm=true (default is false) the result is
    * multiplied by the normalization constant, i.e. divided by
    * sqrt(2*Pi)*sigma.
    */
   double LogGaus(double x, double mean = 0, double sigma = 1, bool norm = false);

   /**
    * Calculate the natural logarithm of a poisson distribution.
    */
   double LogPoisson(double x, double par);

   /**
    * Calculates Binomial probability using approximations for
    * factorial calculations if calculation for number greater than 20
    * required using the BCMath::ApproxLogFact function.
    */
   double ApproxBinomial(int n, int k, double p);

   /**
    * Calculates natural logarithm of the Binomial probability using
    * approximations for factorial calculations if calculation for
    * number greater than 20 required using the BCMath::ApproxLogFact
    * function.
    */
   double LogApproxBinomial(int n, int k, double p);

   /**
    * Calculates natural logarithm of the Binomial factor "n over k"
    * using approximations for factorial calculations if calculation
    * for number greater than 20 required using the
    * BCMath::ApproxLogFact function.  Even for large numbers the
    * calculation is performed precisely, if n-k < 5
    */
   double LogBinomFactor(int n, int k);

   /**
    * Calculates natural logarithm of the n-factorial (n!) using
    * Srinivasa Ramanujan approximation
    * log(n!) = n*log(n) - n + log(n*(1.+4.*n*(1.+2.*n)))/6. + log(PI)/2.
    * if n > 20.  If n <= 20 it uses BCMath::LogFact to calculate it
    * exactly.
    */
   double ApproxLogFact(double x);

   /**
    * Calculates natural logarithm of the Binomial factor "n over k".
    */
   double LogNoverK(int n, int k);

   /**
    * Calculates natural logarithm of the n-factorial (n!)
    */
   double LogFact(int n);

   /**
    * Returns the "greater or equal" of two numbers
    */
   inline int Max(int x, int y)
      { return x >= y ? x : y; }

   inline int Max(unsigned int x, unsigned int y)
      { return x >= y ? x : y; }

   inline double Max(double x, double y)
      { return x >= y ? x : y; }

   inline double Max(float x, float y)
      { return x >= y ? x : y; }

   /**
    * Returns the "less or equal" of two numbers
    */
   inline int Min(int x, int y)
      { return x <= y ? x : y; }

   inline int Min(unsigned int x, unsigned int y)
      { return x <= y ? x : y; }

   inline double Min(double x, double y)
      { return x <= y ? x : y; }

   inline double Min(float x, float y)
      { return x <= y ? x : y; }

   /**
    * Returns the nearest integer of a double number.
    */
   int Nint(double x);

   /**
    * Returns the rms of an array.
    */
   double rms(int n, const double * a);

   /**
    * Calculates the logarithm of the non-relativistic Breit-Wigner
    * distribution.
    */
   double LogBreitWignerNonRel(double x, double mean, double Gamma, bool norm = false);
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
    * Get N random numbers distributed according to chi square function
    * with K degrees of freedom
    */
   void RandomChi2(std::vector<double> &randoms, int K);


   /**
    * Calculate the empirical cumulative distribution function for
    * one dimensional data vector. For consistency, the ECDF
    * of value smaller than the minimum observed (underflow bin) is zero, and
    * for larger than maximum (overflow bin) it is one.
    *
    * @param   data  the observations
    * @return  histogram with normalized ECDF
    */
   TH1D* ECDF(const std::vector<double>& data);

   /**
    * Find the longest runs of zeros and ones in
    * the bit stream
    *
    * @param   bitStream input sequence of boolean values
    * @return  runs  first entry the longest zeros run, second entry the longest ones run
    */
   std::vector<int> longestRuns(const std::vector<bool>& bitStream);

   /**
    * Find the longest success/failure run in set of norm. distributed variables.
    *  Success = observation >= expectation.
    * Runs are weighted by the total chi^2 of all elements in the run
    *
    * @param   yMeasured the observations
    * @param  yExpected   the expected values
    * @param  sigma the theoretical uncertainties on the expectations
    * @return  runs  first entry the max. weight failure run,
    *            second entry the max. success run */
   std::vector<double> longestRunsChi2(const std::vector<double>& yMeasured,
      const std::vector<double>& yExpected, const std::vector<double>& sigma);

   /**
    * Find the sampling probability that, given n independent Bernoulli
    * trials with success rate = failure rate = 1/2, the longest run of
    * consecutive successes is greater than the longest observed run.
    * Key idea from
    * Burr, E.J. & Cane, G. Longest Run of Consecutive Observations Having a Specified Attribute. Biometrika 48, 461-465 (1961).
    *
    *
    * @param   longestObserved  actual longest run
    * @param   nTrials number of independent trials
    * @return  frequency
    */
   double longestRunFrequency(unsigned int longestObserved, unsigned int nTrials);


   double SplitGaussian(double* x, double* par);


   /** Cache factorials for first \arg \c n integers.
    * The cache is filled upon first call of LogFact(). */
   void CacheFactorial(unsigned int n);

	/**
	 * Compute the R-value according to Gelman-Rubin,
	 * GR1992 : Gelman, A. and Rubin, D.B.
     * Inference from Iterative Simulation Using Multiple Sequences
     * Statistical Science, Vol. 7, No. 4 (Nov. 1992), pp. 457-472
	 *
	 * @param chain_means
	 * @param chain_variances
	 * @param chain_length
	 * @param strict If true, use the algorithm laid forth in the paper,
	 *               else use a relaxed version which generally leads to smaller R-values.
	 * @return R-value
	 */
	double Rvalue(const std::vector<double> & chain_means, const std::vector<double> & chain_variances,
                  const unsigned & chain_length, const bool & strict = true) throw (std::invalid_argument, std::domain_error);
}

// ---------------------------------------------------------

#endif

