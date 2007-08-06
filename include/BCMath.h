/**
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
 * --------------------------------------------------------- 
 *
 *
 * A namespace which encapsulates some mathematical functions 
 * necessary for BAT.
 *
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
	/*
	 * Calculate the natural logarithm of a gaussian function with mean and sigma.
	 * If norm=true (default is false) the result is multiplied by the normalization
	 * constant, i.e. divided by sqrt(2*Pi)*sigma.
	 */
	double LogGaus(double x, double mean = 0, double sigma = 1, bool norm = false);
}; 

// --------------------------------------------------------- 

#endif 
