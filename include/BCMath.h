/*!
 * \namespace BCMath
 * \brief Some useful mathematic functions.
 * \author Daniel Kollar
 * \author Kevin Kr&ouml;ninger
 * \version 1.0
 * \date 08.2008
 * \detail A namespace which encapsulates some mathematical functions
 * necessary for BAT.
 */

/*
 * Copyright (C) 2008, Daniel Kollar and Kevin Kroeninger.
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

// ---------------------------------------------------------

#ifndef __BCMATH__H
#define __BCMATH__H

#include <math.h>

#define BCMATH_NFACT_ALIMIT 20

// ---------------------------------------------------------

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

	inline double Max(double x, double y)
		{ return x >= y ? x : y; }

	/**
	 * Returns the "less or equal" of two numbers
	 */
	inline int Min(int x, int y)
		{ return x <= y ? x : y; }

	inline double Min(double x, double y)
	{ return x <= y ? x : y; }

	/**
	 * Returns the nearest integer of a double number.
	 */
	int Nint(double x);

	/**
	 * Returns the rms of an array.
	 */
	double rms(int n, const double * a);

	/*
	 * Calculates the logarithm of the non-relativistic Breit-Wigner
	 * distribution.
	 */
	double LogBreitWignerNonRel(double x, double mean, double Gamma, bool norm = false);
	double LogBreitWignerRel(double x, double mean, double Gamma);

};

// ---------------------------------------------------------

#endif

