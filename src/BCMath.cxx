/*
 * Copyright (C) 2008, Daniel Kollar and Kevin Kroeninger.
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

// ---------------------------------------------------------

#include "BCMath.h"

#include <math.h>

#include <TMath.h>

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

