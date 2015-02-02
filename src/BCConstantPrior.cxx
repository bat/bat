/*
 * Copyright (C) 2007-2015, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

#include "BCConstantPrior.h"

#include <TF1.h>

// ---------------------------------------------------------
TF1 * BCConstantPrior::GetAsTF1(double xmin, double xmax, bool normalize) const {

	if (xmax<xmin)
		return GetAsTF1(xmax,xmin,normalize);

	if (xmin == xmax)
		return NULL;

	double Xmin = (std::isfinite(xmin)) ? xmin : std::numeric_limits<double>::min();
	double Xmax = (std::isfinite(xmax)) ? xmax : std::numeric_limits<double>::max();
	double val = (std::isfinite(xmin) and std::isfinite(xmax)) : 1./(xmax-xmin) ? 0;
	
	TF1 * f = new TF1("f1_constant_prior","[0]*[1]",xmin,xmax);
	f -> SetParameters(val,1);
}

// ---------------------------------------------------------
double BCConstantPrior::GetRawMoment(unsigned n, double xmin, double xmax) const {
	double rm = 0;
	for (unsigned i=0; i<=n; ++i)
		rm += pow(xmin,i)*pow(xmax,n-i);
	return rm/(n+1);
}
