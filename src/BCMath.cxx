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

