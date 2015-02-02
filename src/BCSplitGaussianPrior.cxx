/*
 * Copyright (C) 2007-2015, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

#include "BCGaussianPrior.h"

#include <cmath>

#include <TF1.h>

// ---------------------------------------------------------
TF1 * BCSplitGaussianPrior::GetAsTF1(double xmin, double xmax, bool normalize) const {
	if (xmax<xmin)
		return GetAsTF1(xmax,xmin,normalize);

	if (xmax==xmin)
		return NULL;

	double Xmin = (std::isfinite(xmin)) ? xmin : std::numeric_limits<double>::min();
	double Xmax = (std::isfinite(xmax)) ? xmax : std::numeric_limits<double>::max();

	double integral = (normalize) ? GetIntegral(xmin,xmax) : 1;

	TF1 * f = new TF1("f1_gaussian_prior","sqrt(2/pi)/([1]+[2]) * exp(-0.5*((x-[0])/((x<=[0])*[1]+(x>[0])*[1]))^2)",Xmin,Xmax);
	f -> SetParameters(fMean,fSigmaLow,fSigmaHigh,integral);
	return f;
}

// ---------------------------------------------------------
double BCSplitGaussianPrior::GetLogPrior(double x) const {
	if (x > fMean)
		return -(x-fMean)*(x-fMean)/fSigmaHigh/fSigmaHigh/2. + log(2/M_PI) - log(fSigmaHigh+fSigmaLow);
	return -(x-fMean)*(x-fMean)/fSigmaLow/fSigmaLow/2. + log(2/M_PI) - log(fSigmaHigh+fSigmaLow);
}

// ---------------------------------------------------------
double BCSplitGaussianPrior::GetRawMoment(unsigned n, double xmin, double xmax) const {
	if (n==0)
		return std::numeric_limits<double>::infinity();

	BCPrior::BCPriorRange r = CheckLimits(xmin,xmax);

	if (r==BCPrior::kEmptyRange)
		return (n==1) ? xmin : 0;

	switch (n) {
		
	case 1: {											// mean
		if (r==kInfiniteRange)
			return fMean;

		double phi_min = (r==BCPrior::kNegativeInfiniteRange) ? 0 : fSigmaLow*fSigmaLow   * TMath::Gaus(xmin,fMean,fSigmaLow,false)/sqrt(2*M_PI); 
		double phi_max = (r==BCPrior::kPositiveInfiniteRange) ? 0 : fSigmaHigh*fSigmaHigh * TMath::Gaus(xmax,fMean,fSigmaHigh,false)/sqrt(2*M_PI); 

		double erf_min = (r==kNegativeInfiniteRange) ? -fSigmaLow  : fSigmaLow  * TMath::Erf((xmin-fMean)/fSigmaLow/sqrt(2));
		double erf_max = (r==kPositiveInfiniteRange) ? +fSigmaHigh : fSigmaHigh * TMath::Erf((xmax-fMean)/fSigmaHigh/sqrt(2));

		double delta_s2 = (fSigmaHigh*fSigmaHigh - fSigmaLow*fSigmaLow) * sqrt(2/M_PI);

		if (xmin > fMean) {
			phi_min = fSigmaHigh*fSigmaHigh * TMath::Gaus(xmin,fMean,fSigmaHigh,false)/sqrt(2*M_PI);
			erf_min = fSigmaHigh * TMath::Erf((xmin-fMean)/fSigmaHigh/sqrt(2));
			delta_s2 = 0;
		}
		if (xmax < fMean) {
			phi_max = fSigmaLow*fSigmaLow   * TMath::Gaus(xmax,fMean,fSigmaLow,false)/sqrt(2*M_PI); 
			erf_max = fSigmaLow  * TMath::Erf((xmax-fMean)/fSigmaLow/sqrt(2));
			delta_s2 = 0;
		}

		return fMean + (delta_s2 - 2*(phi_max-phi_min))/(erf_max-erf_min);
	}

	case 2: {											// second moment
		if (r==kInfiniteRange)
			return fMean*fMean + (fSigmaHigh*fSigmaHigh + fSigmaLow*fSigmaLow + (fSigmaHigh-fSigmaLow)*fMean)/2;

		double phi_min = (r==BCPrior::kNegativeInfiniteRange) ? 0 : fSigmaLow*fSigmaLow   * (fMean+xmin) * TMath::Gaus(xmin,fMean,fSigmaLow,false) * sqrt(2/M_PI); 
		double phi_max = (r==BCPrior::kPositiveInfiniteRange) ? 0 : fSigmaHigh*fSigmaHigh * (fMean+xmax) * TMath::Gaus(xmax,fMean,fSigmaHigh,false) *sqrt(2/M_PI); 

		double erf_min = (r==kNegativeInfiniteRange) ? -fSigmaLow  : fSigmaLow  * TMath::Erf((xmin-fMean)/fSigmaLow/sqrt(2));
		double erf_max = (r==kPositiveInfiniteRange) ? +fSigmaHigh : fSigmaHigh * TMath::Erf((xmax-fMean)/fSigmaHigh/sqrt(2));

		double k_min = erf_min*fSigmaLow*fSigmaLow;
		double k_max = erf_max*fSigmaHigh*fSigmaHigh;

		double delta_s2 = 2*fMean*(fSigmaHigh*fSigmaHigh-fSigmaLow*fSigmaLow)*sqrt(2/M_PI);

		if (xmin > fMean) {
			phi_min = fSigmaHigh*fSigmaHigh * (fMean+xmin) * TMath::Gaus(xmin,fMean,fSigmaHigh,false) *sqrt(2/M_PI);
			erf_min = fSigmaHigh * TMath::Erf((xmin-fMean)/fSigmaHigh/sqrt(2));
			k_min = erf_min * fSigmaHigh*fSigmaHigh;
			delta_s2 = 0;
		}
		if (xmax < fMean) {
			phi_max = fSigmaLow*fSigmaLow   * (fMean+xmax) * TMath::Gaus(xmax,fMean,fSigmaLow,false) * sqrt(2/M_PI); 
			erf_max = fSigmaLow  * TMath::Erf((xmax-fMean)/fSigmaLow/sqrt(2));
			k_max = erf_max * fSigmaLow*fSigmaLow;
			delta_s2 = 0;
		}

		return fMean*fMean + ( (k_max-k_min) - delta_s2 - (phi_max-phi_min) ) / (erf_max-erf_min);
	}

	default:
		return std::numeric_limits<double>::infinity();

	}
}

// ---------------------------------------------------------
double BCSplitGaussianPrior::GetIntegral(double xmin, double xmax) const {
	BCPrior::BCPriorRange r = CheckLimits(xmin,xmax);

	if (r==kEmptyRange)
		return 0;

	if (r==kInfiniteRange)
		return 1;

	double erf_min = (r==kNegativeInfiniteRange) ? -fSigmaLow  : fSigmaLow  * TMath::Erf((xmin-fMean)/fSigmaLow/sqrt(2));
	double erf_max = (r==kPositiveInfiniteRange) ? +fSigmaHigh : fSigmaHigh * TMath::Erf((xmax-fMean)/fSigmaHigh/sqrt(2));
	
	if (xmin > fMean)
		erf_min = fSigmaHigh * TMath::Erf((xmin-fMean)/fSigmaHigh/sqrt(2));
	if (xmax < fMean)
		erf_max = fSigmaLow  * TMath::Erf((xmax-fMean)/fSigmaLow/sqrt(2));

	return (erf_max-erf_min)/(fSigmaHigh+fSigmaLow);
}

