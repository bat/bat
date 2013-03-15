// ***************************************************************
// This file was created using the CreateProject.sh script.
// CreateProject.sh is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://www.mppmu.mpg.de/bat
// ***************************************************************

#include "ReferenceCounting.h"

#include <BAT/BCMath.h>

#include <TMath.h>
#include <TH1D.h>

#include <iostream>
#include <iomanip>

// ---------------------------------------------------------
ReferenceCounting::ReferenceCounting() : BCModel()
																			 , fEvalOption(kHistogram)
																			 , fNObs(0)
																			 , fEps(1e-6)
																			 , fAlpha(10)
																			 , fBeta(5)
																			 , logs(0)
																			 , maxn(0)
																			 , maxk(10000)
																			 , fHistPriorS(0)
																			 , fHistPriorB(0)
{
   // default constructor
   DefineParameters();
}

// ---------------------------------------------------------
ReferenceCounting::ReferenceCounting(const char * name) : BCModel(name)
																												, fEvalOption(kHistogram)
																												, fNObs(0)
																												, fEps(1e-6)
																												, fAlpha(10)
																												, fBeta(5)
																												, logs(0)
																												, maxn(0)
																												, maxk(10000)
																												, fHistPriorS(0)
																												, fHistPriorB(0)
{
   // constructor
   DefineParameters();
}

// ---------------------------------------------------------
ReferenceCounting::~ReferenceCounting()
   // default destructor
{
	if (fHistPriorS)
		delete fHistPriorS;

	if (fHistPriorB)
		delete fHistPriorB; 
}

// ---------------------------------------------------------
void ReferenceCounting::SetAlphaBeta(double alpha, double beta)
{
	fAlpha = alpha; 
	fBeta  = beta;

	// clear vectors
	helper_a.clear();
	helper_b.clear();
	helper_c.clear();
}

// ---------------------------------------------------------
void ReferenceCounting::DefineParameters()
{
   // Add parameters to your model here.
   // You can then use them in the methods below by calling the
   // parameters.at(i) or parameters[i], where i is the index
   // of the parameter. The indices increase from 0 according to the
   // order of adding the parameters.

	AddParameter("s", 0, 50.); // index 0
	AddParameter("b", 0, 50.); // index 1
}

// ---------------------------------------------------------
double ReferenceCounting::LogLikelihood(const std::vector<double> &parameters)
{
   // This methods returns the logarithm of the conditional probability
   // p(data|parameters). This is where you have to define your model.

	if (!fHistPriorS)
		FillPriorS();

	if (!fHistPriorB)
		FillPriorB();

   double logprob = 0.;

   double s = parameters.at(0);
   double b = parameters.at(1);
	 double nu = s+b;

	 logprob += BCMath::LogPoisson( double(fNObs), nu);

   return logprob;
}

// ---------------------------------------------------------
double ReferenceCounting::LogAPrioriProbability(const std::vector<double> &parameters)
{
   // This method returns the logarithm of the prior probability for the
   // parameters p(parameters).

   double logprob = 0.;

   double s = parameters.at(0);
   double b = parameters.at(1);

	 if (fEvalOption == kAnalytic) {
		 logprob += log(RefPriorS(s));
		 logprob += log(ConjPriorPoisson(b, fAlpha, fBeta));
	 }
	 else if (fEvalOption == kHistogram) {
		 int bins = fHistPriorS->FindBin(s);
		 int binb = fHistPriorB->FindBin(b);
		 logprob += log( fHistPriorS->GetBinContent(bins) );
		 logprob += log( fHistPriorB->GetBinContent(binb) );
	 }
	 else if (fEvalOption == kHistogram) {
		 // debugKK: Needs to be implemented
		 // ...
	 }

   return logprob;
}

// ---------------------------------------------------------
double ReferenceCounting::ConjPriorPoisson(double nu, double alpha, double beta)
{
	return TMath::GammaDist(nu, alpha, 0., 1/beta);
}

// ---------------------------------------------------------
double ReferenceCounting::RefPriorS(double s)
{
	return sqrtI(s) / sqrtI(0.);
}

// ---------------------------------------------------------
double ReferenceCounting::sqrtI(double s)
{
	double sum = 0;

	double inc = 1.;

	for (int n=0; n<=maxk; ++n) {
		double t1 = f(s, n);
		double t2 = f(s, n+1);
		double dsum = t1*t1/t2;
		sum += dsum;
		inc = fabs( dsum/sum );

		// check if sum does not change by more than predefined precision
		if (inc < fEps)
			break;
	}

	return sqrt( fabs( TMath::Power(fBeta/(1+fBeta), fAlpha) * exp(-s) * sum - 1 ) );
}

// ---------------------------------------------------------
double ReferenceCounting::f(double s, int k)
{
	if (s > 0)
		logs = log(s);
	else {
		if (k==0) 
			return 1;
		else 
			return (fAlpha+k-1)/(k*(1+fBeta)) * f(0., k-1);
	}

	double ftemp = 0;

	extend_abc(k+2);

	for (int n = 0; n <= k; ++n) {
		ftemp += df(s, k, n);
	}

	return ftemp;
}

// ---------------------------------------------------------
double ReferenceCounting::logdf(double s, int k, int n)
{
	double logdf = 0;

	logdf = helper_a[n] + (k-n) * logs - helper_b[n] - helper_c[k-n];

	return logdf;
}

// ---------------------------------------------------------
void ReferenceCounting::extend_abc(int n)
{
	if (n <= maxn)
		return;
	
	int nfrom = maxn;
	
	if (nfrom == 0) {
		helper_a.push_back( log(1) );
		helper_b.push_back( 0. );
    helper_c.push_back( 0. );
		++nfrom;
	}
	if (nfrom == 1) {
		helper_a.push_back( log(fAlpha) );
		helper_b.push_back( log( 1. + fBeta) );
		helper_c.push_back( helper_c[0] );
		++nfrom;
	}
	
	if (maxn > 1)
		nfrom++;

	for (int i = nfrom; i <= n; ++i) {
		double temp_a = helper_a[i-1] + log( 1. + (fAlpha-1.)/double(i)); 
		helper_a.push_back(temp_a);
		helper_b.push_back( int(i) * helper_b[1] );
		helper_c.push_back( log( double(i) ) + helper_c[i-1]); 
	}

	maxn = int( helper_a.size() - 1);
}

// ---------------------------------------------------------
void ReferenceCounting::FillPriorS()
{
	// remove old prior histogram
	if (fHistPriorS)
		delete fHistPriorS;

	// create new histogram
	fHistPriorS = new TH1D("hist_prior_s", ";s;p(s)", fMCMCH1NBins[0], fMCMCBoundaryMin[0], fMCMCBoundaryMax[0]);
	
	// fill histogram
	for (int i = 1; i <= fMCMCH1NBins[0]; ++i){
		double s = fHistPriorS->GetBinCenter(i);
		double p = RefPriorS(s);
		fHistPriorS->SetBinContent(i, p);
	}
}

// ---------------------------------------------------------
void ReferenceCounting::FillPriorB()
{
	// remove old prior histogram
	if (fHistPriorB)
		delete fHistPriorB;

	// create new histogram
	fHistPriorB = new TH1D("hist_prior_b", ";b;p(b)", fMCMCH1NBins[1], fMCMCBoundaryMin[1], fMCMCBoundaryMax[1]);
	
	// fill histogram
	for (int i = 1; i <= fMCMCH1NBins[1]; ++i){
		double b = fHistPriorB->GetBinCenter(i);
		double p = ConjPriorPoisson(b, fAlpha, fBeta);
		fHistPriorB->SetBinContent(i, p);
	}
}

// ---------------------------------------------------------