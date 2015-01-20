/*
 * Copyright (C) 2007-2013, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

// ---------------------------------------------------------
#include "BCParameter.h"
#include "BCLog.h"

#include <string>
#include <limits>

#include <TString.h>
#include <TRandom.h>
#include <TMath.h>


// ---------------------------------------------------------

BCParameter::BCParameter()
	:	BCVariable()
	, fFixed(false)
	, fFixedValue(std::numeric_limits<double>::infinity())
	, fPriorType(BCParameter::kPriorUnset)
	, fPriorParameters(0,0)
	, fPriorContainer(0)
	, fInterpolatePrior(false)
{
	fPrefix = "Parameter";
}

// ---------------------------------------------------------

BCParameter::BCParameter(const BCParameter & other)
	: BCVariable(other)
	, fFixed(other.fFixed)
	, fFixedValue(other.fFixedValue)
	, fPriorType(BCParameter::kPriorUnset)
	, fPriorParameters(0,0)
	, fPriorContainer(0)
	, fInterpolatePrior(false)
{
	CopyPrior(other);
}

// ---------------------------------------------------------

BCParameter::BCParameter(const char * name, double lowerlimit, double upperlimit, const char * latexname)
	: BCVariable(name,lowerlimit,upperlimit,latexname)
	, fFixed(false)
	, fFixedValue(std::numeric_limits<double>::infinity())
	, fPriorContainer(0)
	, fInterpolatePrior(false)
{
	fPrefix = "Parameter";
}

// ---------------------------------------------------------

BCParameter::~BCParameter() {
	if (fPriorContainer)
		delete fPriorContainer;
}

// ---------------------------------------------------------

bool BCParameter::SetPriorConstant() {
	fPriorType = kPriorConstant;

	if (GetRangeWidth() == 0) {
		fPriorParameters.assign(1,std::numeric_limits<double>::infinity());
		return false;
	}

	if (!std::isfinite(GetRangeWidth())) {
		fPriorParameters.assign(1,1);
		return true;
	}

	fPriorParameters.assign(1,-log(GetRangeWidth()));
	return true;
}

// ---------------------------------------------------------

bool BCParameter::SetPrior(const TF1 * const f) {
	// remove existing prior container
	delete fPriorContainer;
	fPriorContainer = 0;
	fPriorType = kPriorUnset;

	if (!f)
		return false;

	fPriorContainer = new TF1(*f);
	fPriorType = kPriorTF1;
	return true;
}

// ---------------------------------------------------------
bool BCParameter::SetPrior(const TH1 * const h, bool interpolate) {
	// remove existing container
	delete fPriorContainer;
	fPriorContainer = 0;
	fPriorType = kPriorUnset;

	if (!h)
		return false;
	
	// check if histogram is 1d
	if (h->GetDimension() != 1) {
		BCLog::OutError(Form("BCParameter::SetPrior : Histogram given forparameter %s is not 1D.",GetName().data()));
		return false;
	}

	fPriorContainer = (TH1*) h->Clone();
	fPriorType = kPriorTH1;

	// normalize:
	double integral = ((TH1*)fPriorContainer) -> Integral("width");
	if (integral <= 0)
		BCLog::OutWarning(Form("BCParameter::SetPrior : Histogram given for parameter %s integrates to zero (or negative).",GetName().data()));
	((TH1*)fPriorContainer) -> Scale(1./integral);

	fInterpolatePrior = interpolate;
	
	return true;
}


// ---------------------------------------------------------
bool BCParameter::CopyPrior(const BCParameter & other) {
	switch (other.fPriorType) {
		
	case kPriorConstant:
		SetPriorConstant();
		return true;

	case kPriorTF1:
		return (other.GetPriorTF1()) ? SetPrior(other.GetPriorTF1()) : false;

	case kPriorTH1:
		return (other.GetPriorTH1()) ? SetPrior(other.GetPriorTH1(), other.fInterpolatePrior) : false;

	case kPriorGaussian:
		return (other.fPriorParameters.size()>=2) ? SetPriorGauss(other.fPriorParameters[0],other.fPriorParameters[1]) : false;

	case kPriorSplitGaussian:
		return (other.fPriorParameters.size()>=3) ? SetPriorGauss(other.fPriorParameters[0],other.fPriorParameters[1],other.fPriorParameters[2]) : false;

	case kPriorUnset:
	default:
		fPriorType = kPriorUnset;
		return true;
	}
}

// ---------------------------------------------------------
bool BCParameter::SetPriorGauss(double mean, double sigma) {
	if (sigma < 0)
		return SetPriorGauss(mean,-sigma);

	if (sigma == 0) {
		BCLog::OutWarning(TString::Format("BCParameter::SetPriorGauss : attempting to set zero-width Gaussian prior; leaving prior unset, but fixing parameter (%s).",GetName().data()));
		fPriorType = kPriorUnset;
		Fix(mean);
		return false;
 	}
	
	if (!std::isfinite(sigma)) {
		BCLog::OutWarning(TString::Format("BCParameter::SetPriorGauss : attempting to set infinite-width Guassian prior; setting prior to constant (%s)",GetName().data()));
		SetPriorConstant();
		return false;
	}

	fPriorType = kPriorGaussian;
	fPriorParameters.assign(2,mean);
	fPriorParameters[1] = sigma;
	return true;
}

// ---------------------------------------------------------
bool BCParameter::SetPriorGauss(double mean, double sigma_below, double sigma_above) {
	if (sigma_below < 0)
		return SetPriorGauss(mean,-sigma_below,sigma_above);
	if (sigma_above < 0)
		return SetPriorGauss(mean,sigma_below,-sigma_above);

	if (sigma_below == 0 and sigma_above == 0) {
		BCLog::OutWarning(TString::Format("BCParameter::SetPriorGauss : attempting to set 0 width Gaussian prior; leaving prior unset, but fixing parameter (%s).",GetName().data()));
		fPriorType = kPriorUnset;
		Fix(mean);
		return false;
	}

	if (!std::isfinite(sigma_below) and !std::isfinite(sigma_above)) {
		BCLog::OutWarning(TString::Format("BCParameter::SetPriorGauss : attempting to set infinite-width Guassian prior; setting prior to constant (%s)",GetName().data()));
		SetPriorConstant();
		return false;
	}
	
	fPriorType = kPriorSplitGaussian;
	fPriorParameters.assign(3,mean);
	fPriorParameters[1] = sigma_below;
	fPriorParameters[2] = sigma_above;
	return true;
	if (fPriorContainer)
		delete fPriorContainer;
}

// ---------------------------------------------------------
 double BCParameter::GetLogPrior(double x) const {
	 if ( fFixed )
		return 0;
	switch (fPriorType) {

	case kPriorConstant:
		return fPriorParameters[0];

	case kPriorTF1:
		if (!GetPriorTF1())
			break;
		return log(GetPriorTF1()->Eval(x));

	case kPriorTH1:
		if (!GetPriorTH1())
			break;
		if (fInterpolatePrior)
			return log(GetPriorTH1()->Interpolate(x));
		return log(GetPriorTH1()->GetBinContent(GetPriorTH1()->FindFixBin(x)));

	case kPriorGaussian:
		return -0.5*(x-fPriorParameters[0])*(x-fPriorParameters[0])/fPriorParameters[1]/fPriorParameters[1] - log(fPriorParameters[1]) - 0.5*log(2*M_PI);

	case kPriorSplitGaussian:
		if (x <= fPriorParameters[0]) {
			if (fPriorParameters[1]==0) // zero below mean
				return -std::numeric_limits<double>::infinity();
			if (!std::isfinite(fPriorParameters[1])) // constant below mean
				return -log(fPriorParameters[0]-fLowerLimit);
			return -0.5*(x-fPriorParameters[0])*(x-fPriorParameters[0])/fPriorParameters[1]/fPriorParameters[1] - log(fPriorParameters[1]) - 0.5*log(2*M_PI);
		}
		if (fPriorParameters[2]==0)	// zero above mean
			return -std::numeric_limits<double>::infinity();
		if (!std::isfinite(fPriorParameters[2])) // constant above mean
			return -log(fUpperLimit-fPriorParameters[0]);
		return -0.5*(x-fPriorParameters[0])*(x-fPriorParameters[0])/fPriorParameters[2]/fPriorParameters[2] - log(fPriorParameters[2]) - 0.5*log(2*M_PI);
		
	case kPriorUnset:
	default:
		break;
	}

	BCLog::OutError(Form("BCParameter::GetLogPrior : no proper prior defined for parameter %s",GetName().data()));
	return -std::numeric_limits<double>::infinity();
 }

// ---------------------------------------------------------
double BCParameter::GetRandomValueAccordingToPrior(TRandom * const R, unsigned N) const {
	if (fFixed)
		return fFixedValue;

	switch(fPriorType) {

	case kPriorConstant: {
		return GetUniformRandomValue(R);
	}

	case kPriorTF1: {
		if (!GetPriorTF1())
			break;
		return GetPriorTF1() -> GetRandom(fLowerLimit,fUpperLimit);
	}

	case kPriorTH1: {
		if (!GetPriorTH1())
			break;
		return GetPriorTH1() -> GetRandom();
	}

	case kPriorGaussian: {
		if (!R)
			break;
		// make at most N attempts to find value in range
		for (unsigned n=0; n<N; ++n) {
			double x = R -> Gaus(fPriorParameters[0],fPriorParameters[1]);
			if (IsWithinLimits(x))
				return x;
		}
		break;
	}
		
	case kPriorSplitGaussian: {
		double integral_below = TMath::Erf((fPriorParameters[0]-fLowerLimit)/fPriorParameters[1]/sqrt(2));
		double integral_above = TMath::Erf((fUpperLimit-fPriorParameters[0])/fPriorParameters[2]/sqrt(2));
		double ratio = integral_below / (integral_below+integral_above);
		// make at most N attempts to find value in range
		for (unsigned n=0; n<N; ++n) {
			double x = 0;
			// select side of mean to be on
			if (R->Rndm() <= ratio) {
				// below mean:
				x = R -> Gaus(fPriorParameters[0],fPriorParameters[1]);
				// force location to below mean
				if (x > fPriorParameters[0])
					x = 2*fPriorParameters[0] - x;
			} else {
				// above mean:
				x = R -> Gaus(fPriorParameters[0],fPriorParameters[2]);
				// force location to above mean
				if (x < fPriorParameters[0])
					x = 2*fPriorParameters[0] - x;
			}
			if (IsWithinLimits(x))
				return x;
		}
		break;
	}

	case kPriorUnset:
	default:
		break;
	}

	BCLog::OutError(TString::Format("BCParameter::GetRandomValueAccordingToPrior : parameter %s has no valid prior set. Returning NaN.",GetName().data()));
	return std::numeric_limits<double>::quiet_NaN();
}

// ---------------------------------------------------------
double BCParameter::GetPriorMean() const {
	switch (fPriorType) {

	case kPriorConstant: {
		if (!std::isfinite(GetRangeWidth()))
			break;
		return GetRangeCenter();
	}

	case kPriorTF1: {
		if (!GetPriorTF1())
			break;
		return GetPriorTF1() -> Mean(fLowerLimit,fUpperLimit);
	}

	case kPriorTH1: {
		if (!GetPriorTH1())
			break;
		return GetPriorTH1() -> GetMean();
	}

	case kPriorGaussian:
	case kPriorSplitGaussian: {
		return fPriorParameters[0];
	}

	case kPriorUnset:
	default:
		break;
	}
	return std::numeric_limits<double>::infinity();
}

// ---------------------------------------------------------
double BCParameter::GetPriorVariance() const {
	switch (fPriorType) {

	case kPriorConstant: {
		if (!std::isfinite(GetRangeWidth()))
			break;
		return GetRangeWidth()*GetRangeWidth()/12.;
	}

	case kPriorTF1: {
		if (!GetPriorTF1())
			break;
		return GetPriorTF1() -> Variance(fLowerLimit,fUpperLimit);
	}

	case kPriorTH1: {
		if (!GetPriorTH1())
			break;
		double s = GetPriorTH1() -> GetRMS();
		return s*s;
	}

	case kPriorGaussian: {
		return fPriorParameters[1]*fPriorParameters[1];
	}

	case kPriorSplitGaussian: {
		double s0 = (std::isfinite(fPriorParameters[1])) ? fPriorParameters[1] : (fPriorParameters[0]-fLowerLimit)/sqrt(12);
		double s1 = (std::isfinite(fPriorParameters[2])) ? fPriorParameters[2] : (fUpperLimit-fPriorParameters[0])/sqrt(12);
		// use maximum variance
		double s = std::max<double>(s0,s1);
		return s*s;
	}		
		
	case kPriorUnset:
	default:
		break;
	}
	return std::numeric_limits<double>::infinity();
}

// ---------------------------------------------------------
double BCParameter::GetPriorStandardDeviation() const {
	switch (fPriorType) {

	case kPriorConstant: {
		if (!std::isfinite(GetRangeWidth()))
			break;
		return GetRangeWidth()/sqrt(12.);
	}

	case kPriorTF1: {
		if (!GetPriorTF1())
			break;
		return sqrt(GetPriorTF1() -> Variance(fLowerLimit,fUpperLimit));
	}

	case kPriorTH1: {
		if (!GetPriorTH1())
			break;
		return GetPriorTH1() -> GetRMS();
	}

	case kPriorGaussian: {
		return fPriorParameters[1];
	}

	case kPriorSplitGaussian: {
		double s0 = (std::isfinite(fPriorParameters[1])) ? fPriorParameters[1] : (fPriorParameters[0]-fLowerLimit)/sqrt(12);
		double s1 = (std::isfinite(fPriorParameters[2])) ? fPriorParameters[2] : (fUpperLimit-fPriorParameters[0])/sqrt(12);
		// use maximum standard deviation
		return std::max<double>(s0,s1);
	}		
		
	case kPriorUnset:
	default:
		break;
	}
	return std::numeric_limits<double>::infinity();
}

// ---------------------------------------------------------
double BCParameter::GetRandomValueAccordingToGaussianOfPrior(TRandom * const R, double expansion_factor, unsigned N) const {
	double m = GetPriorMean();

	// if mean is not finite, return it
	if (!std::isfinite(m))
		return m;
	
	// if no random number generator provided, return nan
	if (!R)
		return std::numeric_limits<double>::quiet_NaN();

	double s = GetPriorStandardDeviation();

	// increase standard deviation by expansion_factor
	s *= expansion_factor;

	// if standard deviation is infinite, use uniform distribution
	if (!std::isfinite(s))
		return GetUniformRandomValue(R);

	// if standard deviation is zero, use delta(x-mean)
	if (s==0)
		return m;

	// check overlap of range with gaussian and issue warning if necessary
	if (IsWithinLimits(m)) {
		if (GetRangeWidth()/s < 1e-4)
			BCLog::OutWarning("BCParameter::GetRandomValueAccordingToGaussianOfPrior : Parameter range very small compared to prior width, this may take a while...");
	} else {
		double d = ((m<fLowerLimit) ? fLowerLimit-m : m-fUpperLimit) / s;
		if (d>4)
			BCLog::OutWarning("BCParameter::GetRandomValueAccordingToGaussianOfPrior : Parameter range beyond 3 sigma of prior mean, this may take a while...");
		else if (GetRangeWidth()/s < 1e-4)
			BCLog::OutWarning("BCParameter::GetRandomValueAccordingToGaussianOfPrior : Parameter range very small compared to prior width, this may take a while...");
	}

	// repeat until value inside range is selected
	// or maximum number of tries is exhausted
	for (unsigned n=0; n<N; ++n) {
		double x = R->Gaus(m,s);
		if (IsWithinLimits(x))
			return x;
	}
	return std::numeric_limits<double>::quiet_NaN();
}

// ---------------------------------------------------------
std::string BCParameter::OneLineSummary() const {
	if (!Fixed())
		return BCVariable::OneLineSummary();
	return std::string(TString::Format("%s (fixed at %.*f)",BCVariable::OneLineSummary().data(),GetPrecision(),GetFixedValue()));
}

