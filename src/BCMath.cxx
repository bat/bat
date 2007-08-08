#include "BCMath.h"

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
	return result - TMath::Log( TMath::Sqrt( 2.*TMath::Pi() ) * sigma );
}

// --------------------------------------------------------- 
double BCMath::ApproxBinomial(int n, int k, double p)
{
	return TMath::Exp( BCMath::LogApproxBinomial(n, k, p) );
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

	return BCMath::LogBinomFactor(n,k) + (double)k*TMath::Log(p) + (double)(n-k)*TMath::Log(1.-p);
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
	if(k==1 || k==n-1) return TMath::Log((double)n);

	// set treshold for using approximations
	int flimit = 100;

	// if no approximation needed
	if(n<flimit || (n-k)<5) return BCMath::LogNoverK(n,k);

	// calculate final log(n over k) using approximations if necessary
	return BCMath::ApproxLogFact((double)n) - BCMath::ApproxLogFact((double)k) - BCMath::ApproxLogFact((double)(n-k));
}

// --------------------------------------------------------- 
double BCMath::ApproxLogFact(double x)
{
	int flimit=100;
	if(x>flimit)
		return x*TMath::Log(x) - x + .5*TMath::Log(2.*TMath::Pi()*x) - 1./(12.*x);
	else
		return BCMath::LogFact((int)x);
}

// --------------------------------------------------------- 
double BCMath::LogFact(int n)
{
	double ln = 0.;
	for(int i=1;i<n;i++)
		ln += TMath::Log((double)i);
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
	if(k==1 || k==n-1) return TMath::Log((double)n);

	int lmax = TMath::Max(k,n-k);
	int lmin = TMath::Min(k,n-k);

	double ln = 0.;
	for(int i=n;i>lmax;i--)
		ln += TMath::Log((double)i);
	ln -= BCMath::LogFact(lmin);

	return ln;
}

