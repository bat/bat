/*
 * Copyright (C) 2008-2011, Daniel Kollar, Kevin Kroeninger and Jing Liu
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
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

#include <Math/PdfFuncMathCore.h>

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
      return LogFact((int) x);
}

// ---------------------------------------------------------
double LogFact(int n)
{
	// return NaN for negative argument
	if (n<0)
		return std::numeric_limits<double>::quiet_NaN();

	// cache the factorials on first call
	if ( !BCMath::logfact && BCMath::nCacheFact>=0 ) {
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

   int lmax = Max(k, n - k);
   int lmin = Min(k, n - k);

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
      if (x + .5 == (double) i && i&1)
         i--;
   }
   else {
      i = int(x - 0.5);
      if (x - 0.5 == double(i) && i&1)
         i++;
   }

   return i;
}

// ---------------------------------------------------------

double rms(int n, const double *a)
{
   if (n <= 0 || !a)
      return 0;

   double sum = 0., sum2 = 0.;

   for (int i = 0; i < n; i++) {
      sum += a[i];
      sum2 += a[i] * a[i];
   }

   double n1 = 1. / (double) n;
   double mean = sum * n1;

   return sqrt(fabs(sum2 * n1 - mean * mean));
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
   int nUnique = uniqueObservations.size();
   double lowerEdges[nUnique];

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
   TH1D * ECDF = new TH1D("ECDF", "Empirical cumulative distribution function", nUnique - 1, lowerEdges);

   // fill the data in to find multiplicities
   for (int i = 0; i < N; ++i)
      ECDF -> Fill(data[i]);

   // just in case, empty the underflow
   ECDF -> SetBinContent(0, 0.0);

   // construct the ecdf
   for (int nBin = 1; nBin <= ECDF->GetNbinsX(); nBin++) {
      double previousBin = ECDF -> GetBinContent(nBin - 1);
      // BCLog::OutDebug(Form("n_%d = %.2f", nBin, ECDF -> GetBinContent(nBin) ));
      // BCLog::OutDebug(Form("previous_%d = %.2f", nBin, previousBin));
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
            maxRunAbove = Max(maxRunAbove, currRun);
         else
            maxRunBelow = Max(maxRunBelow, currRun);
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
      maxRunAbove = Max(maxRunAbove, currRun);
   else
      maxRunBelow = Max(maxRunBelow, currRun);

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
            maxRunAbove = Max(maxRunAbove, currRun);
         else
            maxRunBelow = Max(maxRunBelow, currRun);
         //set flag to run of opposite kind
         aboveRun = !aboveRun;
         //restart at current residual
         currRun = residue * residue;
      }
      //BCLog::OutDebug(Form("maxRunBelow = %g", maxRunBelow));
      //BCLog::OutDebug(Form("maxRunAbove = %g", maxRunAbove));
      //BCLog::OutDebug(Form("currRun = %g", currRun));

   }

   //BCLog::OutDebug(Form("maxRunBelow = %g", maxRunBelow));
   //BCLog::OutDebug(Form("maxRunAbove = %g", maxRunAbove));
   //BCLog::OutDebug(Form("currRun = %g", currRun));

   //check last run
   if (aboveRun)
      maxRunAbove = Max(maxRunAbove, currRun);
   else
      maxRunBelow = Max(maxRunBelow, currRun);

   //BCLog::OutDebug(Form("maxRunBelow = %g", maxRunBelow));
   //BCLog::OutDebug(Form("maxRunAbove = %g", maxRunAbove));

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

   // first method: use the gamma function for the factorials: bit slower and more inaccurate
   // in fact may return NaN for n >= 1000

   //   double Gup, Gdown1, Gdown2, Gdown3;
   //
   //   for (uint r = 0; r <= n; r++) {
   //      conditionalProb = 0.0;
   //      for (uint i = 1; ( i <= n-r+1) && (i <= uint(r / double(Lobs + 1)) ); i++) {
   //
   //         Gup =  TMath::Gamma(1 - i * (Lobs + 1) + n);
   //         Gdown1 = TMath::Gamma(1 + i);
   //         Gdown2 = TMath::Gamma(2 - i + n - r);
   //         Gdown3 = TMath::Gamma(1 - i * (Lobs + 1) + r);
   //
   //         //consider the sign of contribution
   //         Gup = i%2 ? Gup : - Gup;
   //
   //         conditionalProb += Gup/(Gdown1 * Gdown2 * Gdown3);
   //
   ////         printf("G(%d,%d)= %.2f %.2f %.2f %.2f",r,i,Gup, Gdown1, Gdown2, Gdown3);
   //
   ////         prob += TMath::Gamma(1 - i * (Lobs + 1) + n)
   ////               / ( TMath::Gamma(1 + i) * TMath::Gamma(2 - i + n - r)
   ////                     * TMath::Gamma(1 - i * (Lobs + 1) + r) );
   ////         printf("prob inside (i=%d) = %.2f",i, prob);
   //      }
   //      prob += (1 + n - r)*conditionalProb;
   ////      printf("prob outside (r=%d) = %.2f",r, prob);
   //   }

   // alternative using log factorial approximations, is faster and more accurate

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
         //         printf("tempLog inside = %.2f",prob);
      }
      //      printf("tempLog outside = %.2f",prob);
      prob += (1 + n - r) * conditionalProb;
   }

   // Bernoulli probability of each permutation
   prob *= pow(2., -double(n));

   return prob;
}

} // end of namespace BCMath

