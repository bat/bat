/*! \namespace BCMath
 *  \brief Some useful mathematical functions
 *
 * A namespace which encapsulates some mathematical functions
 * necessary for BAT.
 *
 * --------------------------------------------------------- 
 *
 * AUTHOR:  D. Kollar
 *
 * CONTACT: dkollar *at* mppmu *dot* mppmu *dot* de, kroening *at* mppmu *dot* mppmu *dot* de 
 *
 * CREATED: 03.08.2007 by Dano
 * 
 * REVISION: 
 *
 *  08.08.2007  Dano  * added functions for binomial distributions and factorials\n
 *
 * --------------------------------------------------------- 
 */ 

// --------------------------------------------------------- 

#ifndef __BCMATH__H
#define __BCMATH__H

#include <iostream>
#include <fstream> 

#include <math.h>
#include <TMath.h>

// --------------------------------------------------------- 

namespace BCMath
{
	/*!
	 * Calculate the natural logarithm of a gaussian function with mean and sigma.
	 * If norm=true (default is false) the result is multiplied by the normalization
	 * constant, i.e. divided by sqrt(2*Pi)*sigma.
	 */
	double LogGaus(double x, double mean = 0, double sigma = 1, bool norm = false);

	/*!
	 * Calculates Binomial probability using approximations for factorial
	 * calculations if calculation for number greater than 100 required
	 * using the BCMath::ApproxLogFact function.
	 */
	double ApproxBinomial(int n, int k, double p);

	/*!
	 * Calculates natural logarithm of the Binomial probability using
	 * approximations for factorial calculations if calculation for number greater
	 * than 100 required using the BCMath::ApproxLogFact function.
	 */
	double LogApproxBinomial(int n, int k, double p);

	/*!
	 * Calculates natural logarithm of the Binomial factor "n over k" using
	 * approximations for factorial calculations if calculation for number
	 * greater than 100 required using the BCMath::ApproxLogFact function.
	 * Even for large numbers the calculation is performed precisely, if
	 * n-k < 5
	 */
	double LogBinomFactor(int n, int k);

	/*!
	 * Calculates natural logarithm of the n-factorial (n!) using approximation
	 * log(n!) = n*log(n) - n + .5*log(2*pi*n) -1/(12*x) if n > 100.
	 * If n<=100 it uses BCMath::LogFact to calculate it exactly.
	 */
	double ApproxLogFact(double x);

	/*!
	 * Calculates natural logarithm of the Binomial factor "n over k".
	 */
	double LogNoverK(int n, int k);

	/*!
	 * Calculates natural logarithm of the n-factorial (n!)
	 */
	double LogFact(int n);

}; 

// --------------------------------------------------------- 

#endif

