/*
 * Copyright (C) 2008, Daniel Kollar, Kevin Kroeninger and Jing Liu
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

// ---------------------------------------------------------

#include "BAT/BCMath.h"
#include "BAT/BCLog.h"

#include <math.h>

#include <set>

#include <TMath.h>
#include <TF1.h>
#include <TH1D.h>

#include <Math/PdfFuncMathCore.h>

// ---------------------------------------------------------

double BCMath::LogGaus(double x, double mean, double sigma, bool norm)
{
	// if we have a delta function, return fixed value
	if(sigma==0.)
		return 0;

	// if sigma is negative use absolute value
	if(sigma<0.)
		sigma *= -1.;

	double arg = (x-mean)/sigma;
	double result = -.5 * arg * arg;

	// check if we should add the normalization constant
	if(!norm)
		return result;

	// subtract the log of the denominator of the normalization constant
	// and return
	return result - log( sqrt( 2.* M_PI ) * sigma );
}

// ---------------------------------------------------------

double BCMath::LogPoisson(double x, double par)
{
	if (par > 899)
		return BCMath::LogGaus(x,par,sqrt(par),true);

	if (x<0)
		return 0;

	if (x == 0.)
		return  -par;

	return x*log(par)-par-BCMath::ApproxLogFact(x);
}

// ---------------------------------------------------------

double BCMath::ApproxBinomial(int n, int k, double p)
{
	return exp( BCMath::LogApproxBinomial(n, k, p) );
}

// ---------------------------------------------------------

double BCMath::LogApproxBinomial(int n, int k, double p)
{
	// check p
	if (p == 0)
		return -1e99;

	else if (p == 1)
		return 0;

	// switch parameters if n < k
	if(n<k)
	{
		int a=n;
		n=k;
		k=a;
	}

	return BCMath::LogBinomFactor(n,k) + (double)k*log(p) + (double)(n-k)*log(1.-p);
}

// ---------------------------------------------------------

double BCMath::LogBinomFactor(int n, int k)
{
	// switch parameters if n < k
	if(n<k)
	{
		int a=n;
		n=k;
		k=a;
	}

	if(n==k || k==0) return 0.;
	if(k==1 || k==n-1) return log((double)n);

	// if no approximation needed
	if(n<BCMATH_NFACT_ALIMIT || (n-k)<5) return BCMath::LogNoverK(n,k);

	// calculate final log(n over k) using approximations if necessary
	return BCMath::ApproxLogFact((double)n) - BCMath::ApproxLogFact((double)k) - BCMath::ApproxLogFact((double)(n-k));
}

// ---------------------------------------------------------

double BCMath::ApproxLogFact(double x)
{
	if(x>BCMATH_NFACT_ALIMIT)
		return x*log(x) - x + log(x*(1.+4.*x*(1.+2.*x)))/6. + log(M_PI)/2.;

	else
		return BCMath::LogFact((int)x);
}

// ---------------------------------------------------------
double BCMath::LogFact(int n)
{
	double ln = 0.;

	for(int i=1;i<=n;i++)
		ln += log((double)i);

	return ln;
}

// ---------------------------------------------------------

double BCMath::LogNoverK(int n, int k)
{
	// switch parameters if n < k
	if(n<k)
	{
		int a = n;
		n = k;
		k = a;
	}

	if(n==k || k==0) return 0.;
	if(k==1 || k==n-1) return log((double)n);

	int lmax = BCMath::Max(k,n-k);
	int lmin = BCMath::Min(k,n-k);

	double ln = 0.;

	for(int i=n;i>lmax;i--)
		ln += log((double)i);
	ln -= BCMath::LogFact(lmin);

	return ln;
}

// ---------------------------------------------------------

int BCMath::Nint(double x)
{
	// round to integer
	int i;

	if (x >= 0)
	{
		i = (int)(x + .5);

		if (x + .5 == (double)i && i & 1)
			i--;
	}

	else
	{
		i = int(x - 0.5);

		if (x - 0.5 == double(i) && i & 1)
			i++;
	}

	return i;
}

// ---------------------------------------------------------

double BCMath::rms(int n, const double *a)
{
	if (n <= 0 || !a)
		return 0;

	double sum = 0., sum2 = 0.;

	for (int i = 0; i < n; i++)
	{
		sum  += a[i];
		sum2 += a[i] * a[i];
	}

	double n1 = 1./(double)n;
	double mean = sum * n1;

	return sqrt(fabs(sum2 * n1 - mean * mean));
}

// ---------------------------------------------------------

double BCMath::LogBreitWignerNonRel(double x, double mean, double Gamma, bool norm)
{
	double bw = log(Gamma) - log( (x - mean) * (x - mean) + Gamma * Gamma / 4.0);

	if (norm)
		bw -= log(2. * M_PI);

	return bw;
}

// ---------------------------------------------------------

double BCMath::LogBreitWignerRel(double x, double mean, double Gamma)
{
	return -log( (x*x - mean*mean) * (x*x - mean*mean) + mean*mean*Gamma*Gamma );
}

// ---------------------------------------------------------

double BCMath::LogChi2(double x, int n)
{
	if (x<0) {
		BCLog::OutWarning("BCMath::LogChi2 : parameter cannot be negative!");
		return -1e99;
	}

	if (x==0 && n==1) {
		BCLog::OutWarning("BCMath::LogChi2 : returned value is infinity!");
		return 1e99;
	}

	double nOver2 = ((double) n)/2.;

	return (nOver2-1.)*log(x) - x/2. - nOver2*log(2) - log(TMath::Gamma(nOver2));
}

// ---------------------------------------------------------

double BCMath::LogVoigtian(double x, double sigma, double gamma)
{
	if (sigma<=0 || gamma<=0) {
		BCLog::OutWarning("BCMath::LogVoigtian : widths are negative or zero!");
		return -1e99;
	}

	return log(TMath::Voigt(x,sigma,gamma));
}

// ---------------------------------------------------------

//wrapper with signature to construct a TF1
double chi2(double *x, double *par) {
	return ROOT::Math::chisquared_pdf(x[0], par[0]);
}

void BCMath::RandomChi2(std::vector<double> &randoms, int K){

	//fixed upper cutoff to 1000, might be too small
	TF1 *f = new TF1("chi2", chi2, 0.0, 1000, 1);
	f->SetParameter(0, K);
	f->SetNpx(500);
	//uses inverse-transform method
	//fortunately CDF only built once
	for (unsigned int i = 0; i < randoms.size(); i++)
		randoms.at(i) = f->GetRandom();
	delete f;
}

TH1D* BCMath::ECDF(const std::vector<double>& data)
{

   int N = data.size();

   std::set<double> uniqueObservations;

   //sort and filter out multiple instances
   for (int i = 0; i < N; ++i) {
      uniqueObservations.insert(data[i]);
   }

   //extract lower edges for CDF histogram
   int nUnique = uniqueObservations.size();
   double lowerEdges[nUnique];

   //traverse the set
   std::set<double>::iterator iter;
   int counter = 0;
   for (iter = uniqueObservations.begin(); iter != uniqueObservations.end(); iter++) {
      lowerEdges[counter] = *iter;
      counter++;
   }

   //create histogram where
   // lower edge of first bin = min. data
   // upper edge of last bin = max. data
   TH1D* ECDF = new TH1D("ECDF", "Empirical cumulative distribution function",
         nUnique - 1, lowerEdges);

   //fill the data in to find multiplicities
   for (int i = 0; i < N-1; ++i) {
      ECDF -> Fill(data[i]);
   }

   //just in case, empty the underflow
   ECDF -> SetBinContent(0, 0.0);

   //construct the ecdf
   for (int nBin = 1; nBin <= ECDF->GetNbinsX(); nBin++) {
      double previousBin = ECDF -> GetBinContent(nBin - 1);
      BCLog::OutDebug(Form("n_%d = %.2f", nBin, ECDF -> GetBinContent(nBin) ));
      BCLog::OutDebug(Form("previous_%d = %.2f", nBin, previousBin));
      double thisBin = ECDF -> GetBinContent(nBin) / double(N);
      ECDF -> SetBinContent(nBin, thisBin + previousBin);

      ECDF -> SetBinError(nBin, 0.0);
   }

   //set the endpoint to 1, so all larger values are at CDF=1
   ECDF -> SetBinContent( ECDF->GetNbinsX()+1, 1.0);

   return ECDF;

}

