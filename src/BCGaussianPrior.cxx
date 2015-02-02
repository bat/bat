/*
 * Copyright (C) 2007-2015, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

#include "BCGaussianPrior.h"

#include <TF1.h>

// ---------------------------------------------------------
BCGaussianPrior::BCGaussianPrior(double mean, double sigma)
	: BCPrior()
	, fMean(mean)
	, fSigma(sigma)
{}

// ---------------------------------------------------------
BCGaussianPrior::BCGaussianPrior(const BCGaussianPrior & other)
	: BCPrior(other)
	, fMean(other.fMean)
	, fSigma(other.fSigma)
{}

// ---------------------------------------------------------
TF1 * BCGaussianPrior::GetAsTF1(double xmin, double xmax, bool normalize) const {
	if (xmax<xmin)
		return GetAsTF1(xmax,xmin,normalize);

	if (xmax==xmin)
		return NULL;

	double Xmin = (std::isfinite(xmin)) ? xmin : std::numeric_limits<double>::min();
	double Xmax = (std::isfinite(xmax)) ? xmax : std::numeric_limits<double>::max();

	TF1 * f = new TF1("f1_gaussian_prior","(1/(2*pi)^0.5/[1] * exp(((x-[0])/[1])^2)) / [3]",xmin,xmax);
	double integral = 1;
	if (normalize) {
		double erf_min = (std::isfinite(xmin)) ? TMath::Erf((Xmin-fMean)/fSigma/sqrt(2)) : -1;
		double erf_max = (std::isfinite(xmax)) ? TMath::Erf((Xmax-fMean)/fSigma/sqrt(2)) : +1;
		double integral = 0.5*(xmax-xmin);
	}
	
	f -> SetParameters(fMean,fSigma,integral);
	return f;
}

// ---------------------------------------------------------
double BCGaussianPrior::GetRawMoment(unsigned n, double xmin, double xmax) const {
	if (n==0)
		return std::numeric_limits<double>::infinity();

	BCPrior::BCPriorRange r = CheckLimits(xmin,xmax);

	if (r==BCPrior::kEmptyRange)
		return (n==1) ? xmin : 0;

	switch (n) {
		
	case 1: {											// mean
		if (r==kInfiniteRange)
			return fMean;
		double gmin = (r==BCPrior::kNegativeInfiniteRange) ? 0 : TMath::Gaus(xmin,fMean,fSigma,true);
		double gmax = (r==BCPrior::kPositiveInfiniteRange) ? 0 : TMath::Gaus(xmax,fMean,fSigma,true);
		return fMean - fSigma*fSigma*(gmax-gmin)/GetIntegral(xmin,xmax)/2;
	}

	case 2: {											// second moment
		if (r==kInfiniteRange)
			return fMean*fMean + fSigma*fSigma;
		double gmin = (r==BCPrior::kNegativeInfiniteRange) ? 0 : xmin*TMath::Gaus(xmin,fMean,fSigma,true);
		double gmax = (r==BCPrior::kPositiveInfiniteRange) ? 0 : xmax*TMath::Gaus(xmax,fMean,fSigma,true);
		return fMean*GetMean(xmin,xmax) + fSigma*fSigma*(1-(gmax-gmin)/GetIntegral(xmin,xmax));
	}

	default:
		return std::numeric_limits<double>::infinity();

	}
}

// ---------------------------------------------------------
double BCGaussianPrior::GetIntegral(double xmin, double xmax) const {
	switch (CheckLimits(xmin,xmax)) {

	case kFiniteRange:
		return (TMath::Erf((xmax-fMean)/fSigma/sqrt(2)) - TMath::Erf((xmin-fMean)/fSigma/sqrt(2))) / 2;

	case kNegativeInfiniteRange:
		return (1 + TMath::Erf((xmax-fMean)/fSigma/sqrt(2))) / 2;

	case kPositiveInfiniteRange:
		return (1 - TMath::Erf((xmin-fMean)/fSigma/sqrt(2))) / 2;

	case kInfiniteRange:
		return 1;

	case kEmptyRange:
		return 0;

	default:
		return std::numeric_limits<double>::infinity();
	}
}
