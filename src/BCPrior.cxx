/*
 * Copyright (C) 2007-2015, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

// ---------------------------------------------------------

#include "BCPrior.h"

#include <TMath.h>
#include <TH1.h>

// ---------------------------------------------------------
BCPrior::BCPriorRange BCPrior::CheckLimits(double xmin, double xmax) const {
	if (xmin == xmax)
		return BCPrior::kEmptyRange;
	if (std::isfinite(xmin) and std::isfinite(xmax))
		return BCPrior::kFiniteRange;
	if (std::isfinite(xmax))
		return BCPrior::kNegativeInfiniteRange;
	if (std::isfinite(xmin))
		return BCPrior::kPositiveInfiniteRange;
	return BCPrior::kInfiniteRange;
}

// ---------------------------------------------------------
double BCPrior::GetCentralMoment(unsigned n, double xmin, double xmax) const {
	if (n == 0)
		return std::numeric_limits<double>::infinity();
	
	if (n == 1)
		return 0;

	double mean = GetMean(xmin,xmax);
	if (!std::isfinite(mean))
		return std::numeric_limits<double>::infinity();

	double cm = 0;
	for (unsigned i=n; i>1; --i) {
		double rm = GetRawMoment(i,xmin,xmax);
		if (!std::isfinite(rm))
			return std::numeric_limits<double>::infinity();
		cm += TMath::Binomial(n,i) * rm * pow(-mean,n-i);
	}
	cm -= (n-1) * pow(-mean,n);
	return cm;
}

// ---------------------------------------------------------
double BCPrior::GetStandardisedMoment(unsigned n, double xmin, double xmax) const {
	double variance = GetVariance(xmin,xmax);
	if (!std::isfinite(variance))
		return std::numeric_limits<double>::infinity();

	double cm = GetCentralMoment(n,xmin,xmax);
	if (!std::isfinite(cm))
		return std::numeric_limits<double>::infinity();

	return cm / pow(variance,n/2.);
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
		double gmin = (r==BCPrior::kNegativeInfiniteRange) ? 0 : fSigmaLow*fSigmaLow*TMath::Gaus(xmin,fMean,fSigmaLow,true);
		double gmax = (r==BCPrior::kPositiveInfiniteRange) ? 0 : fSigmaHigh*fSigmaHigh*TMath::Gaus(xmax,fMean,fSigmaHigh,true);
		return fMean - (gmax-gmin-(fSigmaHigh-fSigmaLow)/sqrt(2*M_PI))/GetIntegral(xmin,xmax)/2;
	}

	case 2: {											// second moment
		if (r==kInfiniteRange)
			return fMean*fMean + (fSigmaHigh*fSigmaHigh + fSigmaLow*fSigmaLow + (fSigmaHigh-fSigmaLow)*fMean)/2;
		double erf_min = (r==BCPrior::kNegativeInfiniteRange) ? -1 : TMath::Erf((xmin-fMean)/fSigmaLow/sqrt(2))/2;
		double erf_max = (r==BCPrior::kPositiveInfiniteRange) ? +1 : TMath::Erf((xmax-fMean)/fSigmaHigh/sqrt(2))/2;
		double gmin = fSigmaLow*(fSigmaLow*erf_min - (((r==BCPrior::kNegativeInfiniteRange) ? 0 : xmin*TMath::Gaus(xmin,fMean,fSigmaLow,false)) - fMean/sqrt(2*M_PI)));
		double gmax = fSigmaHigh*(fSigmaHigh*erf_max - (((r==BCPrior::kPositiveInfiniteRange) ? 0 : xmax*TMath::Gaus(xmax,fMean,fSigmaHigh,false)) - fMean/sqrt(2*M_PI)));
		return fMean*GetMean(xmin,xmax) + (gmax-gmin)/(erf_max-erf_min);
	}

	default:
		return std::numeric_limits<double>::infinity();

	}
}

// ---------------------------------------------------------
double BCSplitGaussianPrior::GetIntegral(double xmin, double xmax) const {
	switch (CheckLimits(xmin,xmax)) {

	case kFiniteRange:
		return (TMath::Erf((xmax-fMean)/fSigmaHigh/sqrt(2)) - TMath::Erf((xmin-fMean)/fSigmaLow/sqrt(2))) / 2;

	case kNegativeInfiniteRange:
		return (1 + TMath::Erf((xmax-fMean)/fSigmaHigh/sqrt(2))) / 2;

	case kPositiveInfiniteRange:
		return (1 - TMath::Erf((xmin-fMean)/fSigmaLow/sqrt(2))) / 2;

	case kInfiniteRange:
		return 1;

	case kEmptyRange:
		return 0;

	default:
		return std::numeric_limits<double>::infinity();
	}
}


// ---------------------------------------------------------
double BCCauchyPrior::GetRawMoment(unsigned n, double xmin, double xmax) const {
	if (n==0)
		return std::numeric_limits<double>::infinity();

	BCPrior::BCPriorRange r = CheckLimits(xmin,xmax);

	if (r==BCPrior::kEmptyRange)
		return (n==1) ? xmin : 0;

	switch (n) {
		
	case 1: {											// mean
		if (r==kInfiniteRange)
			return fMean;
		if (r==kNegativeInfiniteRange)
			return -std::numeric_limits<double>::infinity();
		if (r==kPositiveInfiniteRange)
			return +std::numeric_limits<double>::infinity();
		double L = (xmin-fMean)/fScale;
		double H = (xmax-fMean)/fScale;
		return fMean + fScale/(2*M_PI)*log((1+H*H)/(1+L*L))/GetIntegral(xmin,xmax);
	}

	case 2: {											// second moment
		if (r==kInfiniteRange or r==kNegativeInfiniteRange or r==kPositiveInfiniteRange)
			return std::numeric_limits<double>::infinity();
		return 2*fMean*GetMean(xmin,xmax) - (fMean*fMean+fScale*fScale) + fScale/M_PI*(xmax-xmin)/GetIntegral(xmin,xmax);
	}

	default:
		return std::numeric_limits<double>::infinity();

	}
}

// ---------------------------------------------------------
double BCCauchyPrior::GetIntegral(double xmin, double xmax) const {
	switch (CheckLimits(xmin,xmax)) {

	case kFiniteRange:
		return (atan((xmax-fMean)/fScale)-atan((xmin-fMean)/fScale))/M_PI;

	case kNegativeInfiniteRange:
		return 0.5 + atan((xmax-fMean)/fScale)/M_PI;

	case kPositiveInfiniteRange:
		return 0.5 - atan((xmin-fMean)/fScale)/M_PI;

	case kInfiniteRange:
		return 1;

	case kEmptyRange:
		return 0;

	default:
		return std::numeric_limits<double>::infinity();
	}
}


// ---------------------------------------------------------
double BCTF1Prior::GetRawMoment(unsigned n, double xmin, double xmax) const {
	if (n==0 or !fPriorFunction)
		return std::numeric_limits<double>::infinity();

	return fPriorFunction -> Moment(static_cast<double>(n),xmin,xmax);
}


// ---------------------------------------------------------
BCTH1Prior::BCTH1Prior(TH1 * const h, bool interpolate)
	: BCPrior(BCPrior::kPriorTH1)
	, fPriorHistogram(0)
	, fInterpolate(interpolate)
{
	if (h and h->GetDimension()==1)
		fPriorHistogram = (TH1*) h->Clone();
	double integral = fPriorHistogram -> Integral("width");
	if (integral!=0)
		fPriorHistogram -> Scale(1./integral);
}

// ---------------------------------------------------------
double BCTH1Prior::GetRawMoment(unsigned n, double xmin, double xmax) const {
	if (n==0)
		return 0;
	if (n==1)
		return GetMean(xmin,xmax);
	if (n==2)
		return GetMean(xmin,xmax)*GetMean(xmin,xmax) + GetVariance(xmin,xmax);
	return std::numeric_limits<double>::infinity();
}

// ---------------------------------------------------------
double BCTH1Prior::GetCentralMoment(unsigned n, double xmin, double xmax) const {
	double s = GetStandardDeviation(xmin,xmax);
	if (!std::isfinite(s))
		return std::numeric_limits<double>::infinity();
	return pow(s,n) * GetStandardisedMoment(n,xmin,xmax);
}

// ---------------------------------------------------------
double BCTH1Prior::GetStandardisedMoment(unsigned n, double xmin, double xmax) const {
	if (n==0)
		return 0;
	if (n==1)
		return 1;
	if (n==2)
		return GetSkewness(xmin,xmax);
	if (n==3)
		return GetKurtosis(xmin,xmax);
	return std::numeric_limits<double>::infinity();
}

// ---------------------------------------------------------
double BCTH1Prior::GetMean(double xmin, double xmax) const {
	return (fPriorHistogram) ? GetMean() : std::numeric_limits<double>::infinity();
}

// ---------------------------------------------------------
double BCTH1Prior::GetVariance(double xmin, double xmax) const {
	double s = GetStandardDeviation(xmin,xmax);
	return s*s;
}

// ---------------------------------------------------------
double BCTH1Prior::GetStandardDeviation(double xmin, double xmax) const {
	return (fPriorHistogram) ? fPriorHistogram->GetRMS() : std::numeric_limits<double>::infinity();
}

// ---------------------------------------------------------
double BCTH1Prior::GetSkewness(double xmin, double xmax) const {
	return (fPriorHistogram) ? fPriorHistogram->GetSkewness() : std::numeric_limits<double>::infinity();
}

// ---------------------------------------------------------
double BCTH1Prior::GetKurtosis(double xmin, double xmax) const {
	return (fPriorHistogram) ? fPriorHistogram->GetKurtosis() : std::numeric_limits<double>::infinity();
}
