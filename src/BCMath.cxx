/*
 * Copyright (C) 2007-2014, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

// ---------------------------------------------------------

#include "BCMath.h"
#include "BCLog.h"

#include <math.h>
#include <limits>
#include <set>

#include <TMath.h>
#include <TF1.h>
#include <TH1D.h>
#include <TRandom3.h>

#include <Math/PdfFuncMathCore.h>
#include <Math/QuantFuncMathCore.h>

namespace BCMath {

static unsigned int nCacheFact = 10000;
static double * logfact = 0;

// ---------------------------------------------------------

double LogGaus(double x, double mean, double sigma, bool norm)
{
   // if we have a delta function, return fixed value
   if (sigma == 0.)
      return 0;

   // if sigma is negative use absolute value
   if (sigma < 0.)
      sigma *= -1.;

   double arg = (x - mean) / sigma;
   double result = -.5 * arg * arg;

   // check if we should add the normalization constant
   if (!norm)
      return result;

   // subtract the log of the denominator of the normalization constant
   // and return
   return result - log(sqrt(2. * M_PI) * sigma);
}

// ---------------------------------------------------------

double LogPoisson(double x, double par)
{
   if (par > 899)
      return LogGaus(x, par, sqrt(par), true);

   if (x < 0)
      return 0;

   if (x == 0.)
      return -par;

   return x * log(par) - par - ApproxLogFact(x);
}

// ---------------------------------------------------------

double ApproxBinomial(int n, int k, double p)
{
   return exp(LogApproxBinomial(n, k, p));
}

// ---------------------------------------------------------

double LogApproxBinomial(int n, int k, double p)
{
   // check p
   if (p == 0)
      return -1e99;

   else if (p == 1)
      return 0;

   // switch parameters if n < k
   if (n < k) {
      int a = n;
      n = k;
      k = a;
   }

   return LogBinomFactor(n, k) + (double) k * log(p) + (double) (n - k) * log(1. - p);
}

// ---------------------------------------------------------

double LogBinomFactor(int n, int k)
{
   // switch parameters if n < k
   if (n < k) {
      int a = n;
      n = k;
      k = a;
   }

   if (n == k || k == 0)
      return 0.;
   if (k == 1 || k == n - 1)
      return log((double) n);

   // if no approximation needed
   if ( n < BCMATH_NFACT_ALIMIT || (n < (int) BCMath::nCacheFact &&  (n - k) < 10) )
      return LogNoverK(n, k);

   // calculate final log(n over k) using approximations if necessary
   return ApproxLogFact((double)n) - ApproxLogFact((double)k) - ApproxLogFact((double)(n - k));
}

// ---------------------------------------------------------

double ApproxLogFact(double x)
{
   if (x > BCMath::nCacheFact)
      return x * log(x) - x + log(x * (1. + 4. * x * (1. + 2. * x))) / 6. + log(M_PI) / 2.;

   else
      return LogFact(static_cast<int>(floor(x + 0.5f)));
}

// ---------------------------------------------------------
double LogFact(int n)
{
   // return NaN for negative argument
   if (n<0)
      return std::numeric_limits<double>::quiet_NaN();

   // cache the factorials on first call
   if ( !BCMath::logfact) {
      BCMath::logfact = new double[BCMath::nCacheFact+1];
      double tmplogfact = 0;
      BCMath::logfact[0] = tmplogfact;
      for (unsigned int i=1; i<=BCMath::nCacheFact; i++) {
         tmplogfact += log((double) i);
         BCMath::logfact[i] = tmplogfact;
      }
   }

   // return cached value if available
   if (n <= (int) BCMath::nCacheFact)
      return BCMath::logfact[n];

   // calculate factorial starting from the highest cached value
   double ln(0.);
   if (BCMath::logfact)
      ln = BCMath::logfact[nCacheFact];

   for (int i = BCMath::nCacheFact+1; i <= n; i++)
      ln += log((double) i);

   return ln;
}

// ---------------------------------------------------------

void CacheFactorial(unsigned int n)
{
   nCacheFact = n;
}

// ---------------------------------------------------------

double LogNoverK(int n, int k)
{
   // switch parameters if n < k
   if (n < k) {
      int a = n;
      n = k;
      k = a;
   }

   if (n == k || k == 0)
      return 0.;
   if (k == 1 || k == n - 1)
      return log((double) n);

   int lmax = std::max(k, n - k);
   int lmin = std::min(k, n - k);

   double ln = 0.;

   for (int i = n; i > lmax; i--)
      ln += log((double) i);
   ln -= LogFact(lmin);

   return ln;
}

// ---------------------------------------------------------

int Nint(double x)
{
   // round to integer
   int i;

   if (x >= 0) {
      i = (int) (x + .5);
      if (x + .5 == (double) i && (i&1))
         i--;
   }
   else {
      i = int(x - 0.5);
      if (x - 0.5 == double(i) && (i&1))
         i++;
   }

   return i;
}

// ---------------------------------------------------------

double LogBreitWignerNonRel(double x, double mean, double Gamma, bool norm)
{
   double bw = log(Gamma) - log((x - mean) * (x - mean) + Gamma*Gamma / 4.);

   if (norm)
      bw -= log(2. * M_PI);

   return bw;
}

// ---------------------------------------------------------

double LogBreitWignerRel(double x, double mean, double Gamma)
{
   return -log((x*x - mean*mean) * (x*x - mean*mean) + mean*mean * Gamma*Gamma);
}

// ---------------------------------------------------------

double LogChi2(double x, int n)
{
   if (x < 0) {
      BCLog::OutWarning("BCMath::LogChi2 : parameter cannot be negative!");
      return -1e99;
   }

   if (x == 0 && n == 1) {
      BCLog::OutWarning("BCMath::LogChi2 : returned value is infinity!");
      return 1e99;
   }

   double nOver2 = ((double) n) / 2.;

   return (nOver2 - 1.) * log(x) - x / 2. - nOver2 * log(2) - log(TMath::Gamma(nOver2));
}

// ---------------------------------------------------------
double LogVoigtian(double x, double sigma, double gamma)
{
   if (sigma <= 0 || gamma <= 0) {
      BCLog::OutWarning("BCMath::LogVoigtian : widths are negative or zero!");
      return -1e99;
   }

   return log(TMath::Voigt(x, sigma, gamma));
}

// ---------------------------------------------------------
double LogGammaPDF(double x, double alpha, double beta)
{
	return alpha * log(beta) - TMath::LnGamma(alpha) + (alpha - 1) * log(x) - beta * x;
}

// ---------------------------------------------------------
double LogLogNormal(double x, double mean, double sigma)
{
   // if we have a delta function, return fixed value
   if (sigma == 0.)
      return 0;

   // if sigma is negative use absolute value
   if (sigma < 0.)
      sigma *= -1.;

   double arg = (log(x) - mean) / sigma;
   double result = -.5 * arg * arg;

   return result - log(x * sqrt(2. * M_PI) * sigma);
}

// ---------------------------------------------------------
double SplitGaussian(double* x, double* par)
{
   double mean = par[0];
   double sigmadown = par[1];
   double sigmaup = par[2];

   double sigma = sigmadown;

   if (x[0] > mean)
      sigma = sigmaup;

   return 1.0/sqrt(2.0*M_PI)/sigma * exp(- (x[0]-mean)*(x[0]-mean)/2./sigma/sigma);
}

// ---------------------------------------------------------
// wrapper with signature to construct a TF1
double chi2(double *x, double *par)
{
   return ROOT::Math::chisquared_pdf(x[0], par[0]);
}

// ---------------------------------------------------------
void RandomChi2(std::vector<double> &randoms, int K)
{
   // fixed upper cutoff to 1000, might be too small
   TF1 *f = new TF1("chi2", chi2, 0.0, 1000, 1);
   f->SetParameter(0, K);
   f->SetNpx(500);
   // uses inverse-transform method
   // fortunately CDF only built once
   for (unsigned int i = 0; i < randoms.size(); i++)
      randoms.at(i) = f->GetRandom();

   delete f;
}

// ---------------------------------------------------------
TH1D * ECDF(const std::vector<double> & data)
{
   int N = data.size();

   std::set<double> uniqueObservations;
   // sort and filter out multiple instances
   for (int i = 0; i < N; ++i)
      uniqueObservations.insert(data[i]);

   // extract lower edges for CDF histogram
   unsigned nUnique = uniqueObservations.size();
   std::vector<double> lowerEdges(nUnique);

   // traverse the set
   std::set<double>::iterator iter;
   int counter = 0;
   for (iter = uniqueObservations.begin(); iter != uniqueObservations.end(); ++iter) {
      lowerEdges[counter] = *iter;
      counter++;
   }

   // create histogram where
   // lower edge of first bin = min. data
   // upper edge of last bin = max. data
   TH1D * ECDF = new TH1D("ECDF", "Empirical cumulative distribution function", nUnique - 1, &lowerEdges[0]);

   // fill the data in to find multiplicities
   for (int i = 0; i < N; ++i)
      ECDF -> Fill(data[i]);

   // just in case, empty the underflow
   ECDF -> SetBinContent(0, 0.0);

   // construct the ecdf
   for (int nBin = 1; nBin <= ECDF->GetNbinsX(); nBin++) {
      double previousBin = ECDF -> GetBinContent(nBin - 1);
      double thisBin = ECDF -> GetBinContent(nBin) / double(N);
      ECDF -> SetBinContent(nBin, thisBin + previousBin);

      // the uncertainty is only correctly estimated in the model
      ECDF -> SetBinError(nBin, 0.0);
   }

   // set the endpoint to 1, so all larger values are at CDF=1
   ECDF -> SetBinContent(ECDF->GetNbinsX() + 1, 1.);

   // adjust for nice plotting
   ECDF -> SetMinimum(0.);
   ECDF -> SetMaximum(1.);

   return ECDF;
}

// ---------------------------------------------------------

std::vector<int> longestRuns(const std::vector<bool> &bitStream)
{
   // initialize counter variables
   unsigned int maxRunAbove, maxRunBelow, currRun;
   maxRunAbove = 0;
   maxRunBelow = 0;
   currRun = 1;
   // set both entries to zero
   std::vector<int> runs(2, 0);

   if (bitStream.empty())
      return runs;

   // flag about kind of the currently considered run
   bool aboveRun = bitStream.at(0);

   // start at second variable
   std::vector<bool>::const_iterator iter = bitStream.begin();
   ++iter;
   while (iter != bitStream.end()) {

      // increase counter if run continues
      if (*(iter - 1) == *iter)
         currRun++;
      else {
         // compare terminated run to maximum
         if (aboveRun)
            maxRunAbove = std::max(maxRunAbove, currRun);
         else
            maxRunBelow = std::max(maxRunBelow, currRun);
         // set flag to run of opposite kind
         aboveRun = !aboveRun;
         // restart at length one
         currRun = 1;
      }
      // move to next bit
      ++iter;
   }

   // check last run
   if (aboveRun)
      maxRunAbove = std::max(maxRunAbove, currRun);
   else
      maxRunBelow = std::max(maxRunBelow, currRun);

   // save the longest runs
   runs.at(0) = maxRunBelow;
   runs.at(1) = maxRunAbove;

   return runs;
}
// ---------------------------------------------------------

std::vector<double> longestRunsChi2(
      const std::vector<double>& yMeasured,
      const std::vector<double>& yExpected, const std::vector<double>& sigma)
{
   //initialize counter variables
   double maxRunAbove, maxRunBelow, currRun;
   maxRunAbove = 0;
   maxRunBelow = 0;
   currRun = 0;
   //set both entries to zero
   std::vector<double> runs(2, 0);

   //check input size
   if (yMeasured.size() != yExpected.size() || yMeasured.size() != sigma.size()
         || yExpected.size() != sigma.size()) {
      //should throw exception
      return runs;
   }

   //exclude zero uncertainty
   //...

    int N = yMeasured.size();
    if ( N<=0)
       return runs;
   //BCLog::OutDebug(Form("N = %d", N));


   //flag about kind of the currently considered run
   double residue = (yMeasured.at(0) - yExpected.at(0)) / sigma.at(0);
   bool aboveRun = residue >= 0 ? true : false;
   currRun = residue * residue;

   //start at second variable
   for (int i = 1; i < N; i++) {
      residue = (yMeasured.at(i) - yExpected.at(i)) / sigma.at(i);
      //run continues
      if ((residue >= 0) == aboveRun) {
         currRun += residue * residue;
      } else {
         //compare terminated run to maximum
         if (aboveRun)
            maxRunAbove = std::max(maxRunAbove, currRun);
         else
            maxRunBelow = std::max(maxRunBelow, currRun);
         //set flag to run of opposite kind
         aboveRun = !aboveRun;
         //restart at current residual
         currRun = residue * residue;
      }
   }

   //check last run
   if (aboveRun)
      maxRunAbove = std::max(maxRunAbove, currRun);
   else
      maxRunBelow = std::max(maxRunBelow, currRun);

   //save the longest runs
   runs.at(0) = maxRunBelow;
   runs.at(1) = maxRunAbove;

   return runs;
}
// ---------------------------------------------------------
double longestRunFrequency(unsigned longestObserved, unsigned int nTrials)
{
   // can't observe run that's longer than the whole sequence
   if (longestObserved >= nTrials)
      return 0.;

   // return value
   double prob = 0.;

   // short cuts
   typedef unsigned int uint;
   uint Lobs = longestObserved;
   uint n = nTrials;

   // the result of the inner loop is the cond. P given r successes
   double conditionalProb;

   /* first method: use the gamma function for the factorials: bit slower and more inaccurate
    * in fact may return NaN for n >= 1000.
    * alternative using log factorial approximations, is faster and more accurate
    */

   double tempLog = 0;
   for (uint r = 0; r <= n; r++) {
      conditionalProb = 0.0;

      for (uint i = 1; (i <= n - r + 1) && (i <= uint(r / double(Lobs + 1))); i++) {
         tempLog = ApproxLogFact(n - i * (Lobs + 1)) - ApproxLogFact(i)
               - ApproxLogFact(n - r - i + 1)
               - ApproxLogFact(r - i * (Lobs + 1));
         if (i % 2)
            conditionalProb += exp(tempLog);
         else
            conditionalProb -= exp(tempLog);
      }
      prob += (1 + n - r) * conditionalProb;
   }

   // Bernoulli probability of each permutation
   prob *= pow(2., -double(n));

   return prob;
}

// ---------------------------------------------------------
double Rvalue(const std::vector<double> & chain_means, const std::vector<double> & chain_variances,
              const unsigned & chain_length, const bool & strict)  throw (std::invalid_argument, std::domain_error)
{
   if (chain_means.size() != chain_variances.size())
      throw std::invalid_argument("BCMath::RValue: chain means and chain variances are not aligned!");

   const double n = chain_length;
   const double m = chain_means.size();

   if (m <= 1)
      throw std::invalid_argument("BCMath:RValue: Need at least two chains to compute R-value!");

   // init
   double variance_of_means = 0.0;
   double mean_of_means = 0.0;
   double mean_of_variances = 0.0;
   double variance_of_variances = 0.0;

   // use Welford's method here as well with temporary variables
   double previous_mean_of_means = 0;
   double previous_mean_of_variances = 0;
   for (unsigned i = 0 ; i < m ; ++i)
   {
      if (0 == i)
      {
         mean_of_means = chain_means.front();
         variance_of_means = 0;
         mean_of_variances = chain_variances.front();
         variance_of_variances = 0;

         continue;
      }

      // temporarily store previous mean of step (i-1)
      previous_mean_of_means = mean_of_means;
      previous_mean_of_variances = mean_of_variances;

      // update step
      mean_of_means += (chain_means[i] - previous_mean_of_means) / (i + 1.0);
      variance_of_means += (chain_means[i] - previous_mean_of_means) * (chain_means[i] - mean_of_means);
      mean_of_variances += (chain_variances[i] - previous_mean_of_variances) / (i + 1.0);
      variance_of_variances += (chain_variances[i] - previous_mean_of_variances) * (chain_variances[i] - mean_of_variances);
   }

   variance_of_means /= m - 1.0;
   variance_of_variances /= m - 1.0;

   //use Gelman/Rubin notation
   double B = variance_of_means * n;
   double W = mean_of_variances;
   double sigma_squared = (n - 1.0) / n * W + B / n;

   // avoid case with no variance whatsoever
   if (0.0 == W && 0.0 == B)
   {
      BCLog::OutDebug("BCMath::Rvalue: All samples in all chains identical!");
      return 1.0;
   }
   // avoid NaN due to divide by zero
   if (0.0 == W)
   {
      BCLog::OutDebug("BCMath::Rvalue: W = 0. Avoiding R = NaN.");
      return std::numeric_limits<double>::max();
   }

   // the lax implementation stops here
   if (!strict)
      return sqrt(sigma_squared / W);

   //estimated scale reduction
   double R = 0.0;

   // compute covariances using the means from above
   double covariance_22 = 0.0; // cov(s_i^2, \bar{x_i}^2
   double covariance_21 = 0.0; // cov(s_i^2, \bar{x_i}

   for (unsigned i = 0 ; i < m ; ++i)
   {
      covariance_21 += (chain_variances[i] - mean_of_variances) * (chain_means.at(i) - mean_of_means);
      covariance_22 += (chain_variances[i] - mean_of_variances) * (chain_means[i] * chain_means[i] - mean_of_means * mean_of_means);
   }

   covariance_21 /= m - 1.0;
   covariance_22 /= m - 1.0;

   // scale of t-distribution
   double V = sigma_squared + B / (m * n);

   // estimation of scale variance
   double a = (n - 1.0) * (n - 1.0) / (n * n * m) * variance_of_variances;
   double b = (m + 1) * (m + 1) / (m * n * m * n) * 2.0 / (m - 1) * B * B;
   double c = 2.0 * (m + 1.0) * (n - 1.0) / (m * n * n) * n / m * (covariance_22 - 2.0 * mean_of_means * covariance_21);
   double variance_of_V = a + b + c;

   // degrees of freedom of t-distribution
   double df = 2.0 * V * V / variance_of_V;

   if (df <= 2)
   {
      BCLog::OutDebug(Form("BCMath::Rvalue: DoF (%g) below 2. Avoiding R = NaN.", df));
      return std::numeric_limits<double>::max();;
   }

   // sqrt of estimated scale reduction if sampling were continued
   R = sqrt(V / W * df / (df - 2.0));

   // R smaller, but close to 1 is OK.
   if (R < 0.99 && n > 100)
      throw std::domain_error(Form("BCMath::Rvalue: %g < 0.99. Check for a bug in the implementation!", R));

   return R;
}

// ---------------------------------------------------------
double CorrectPValue(const double & pvalue, const unsigned & npar, const unsigned & nobservations) throw (std::domain_error)
{
   // avoid pathologies
   if (pvalue < 0 or pvalue > 1)
      throw std::domain_error(Form("BCMath::CorrectPValue: pvalue (%g) out of range", pvalue));

   if (pvalue < 1e-150)
      return 0;

   // dof = nobs - npar
   if (npar >= nobservations)
       throw std::domain_error(Form("BCMath::CorrectPValue: "
             "npar exceeds nobservations, %u vs %u", npar, nobservations));

   // convert upper-tail p value to a chi squared
   const double chi2 = ROOT::Math::chisquared_quantile_c(pvalue, nobservations);

   // corrected degrees of freedom
   const unsigned dof = nobservations - npar;

   // transform back to p value
   return TMath::Prob(chi2, dof);
}

// ---------------------------------------------------------
double FastPValue(const std::vector<unsigned> & observed, const std::vector<double> & expected,
                  unsigned nIterations, unsigned seed) throw (std::invalid_argument)
{
   size_t nbins = observed.size();
   if (nbins != expected.size())
   {
      throw std::invalid_argument(Form("BCMath::FastPValue: "
             "size of expected and observed do not match, %u vs %u", unsigned(expected.size()), unsigned(nbins)));
   }

   // temporary histogram to be modified in each MCMC step
   std::vector<unsigned> histogram(nbins, 0);

   // fix seed to iterations for reproducible results
   TRandom3 rng(seed);

   // keep track of log of probability and count data sets with larger value
   double logp = 0;
   double logp_start = 0;
   int counter_pvalue = 0;

   // define starting distribution as histogram with most likely entries
   for (size_t ibin = 0; ibin < nbins; ++ibin) {

      // get the number of expected events
      double yexp = expected[ibin];

      //most likely observed value = int(expected value)
      histogram[ibin] = size_t(yexp);

      // calculate test statistic (= likelihood of entire histogram) for starting distr.
      logp += LogPoisson(size_t(yexp), yexp);

      //statistic of the observed data set
      logp_start += LogPoisson(observed[ibin], yexp);
   }

   // loop over iterations
    for (unsigned iiter = 0; iiter < nIterations; ++iiter) {
       // loop over bins
       for (size_t ibin = 0; ibin < nbins; ++ibin) {
          // random step up or down in statistics for this bin
          double ptest = rng.Rndm() - 0.5;

          // increase statistics by 1
          if (ptest > 0) {
             // calculate factor of probability
             double r = expected[ibin] / double(histogram[ibin] + 1);

             // walk, or don't (this is the Metropolis part)
             if (rng.Rndm() < r) {
                ++histogram[ibin];
                logp += log(r);
             }
          }

          // decrease statistics by 1
          else if (ptest <= 0 && histogram[ibin] > 0) {
             // calculate factor of probability
             double r = double(histogram[ibin]) / expected[ibin];

             // walk, or don't (this is the Metropolis part)
             if (rng.Rndm() < r) {
                --histogram[ibin];
                logp += log(r);
             }
          }
       } // end of looping over bins

       // increase counter
       if (logp <= logp_start)
          ++counter_pvalue;
       // handle the case where a -b +b > a because of precision loss
       else if (logp_start && logp != logp_start && fabs((logp - logp_start) / logp_start) < 1e-15)
          ++counter_pvalue;

    } // end of looping over iterations

      // calculate p-value
      return double(counter_pvalue) / nIterations;
}
} // end of namespace BCMath
