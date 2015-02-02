#ifndef __BCPRIOR__H
#define __BCPRIOR__H

/*!
 * \class BCPrior
 * \brief A class to represent the prior of a parameter
 * \author Daniel Greenwald
 * \version 1.0
 * \date 01.2015
 */


/*
 * Copyright (C) 2007-2015, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

// ---------------------------------------------------------

#include <limits>
#include <cmath>

class TF1;

// ---------------------------------------------------------

class BCPrior {

public:

	/** \name Enumerators */
	/** @{ */

	/** Range types. */
	enum BCPriorRange {
		kFiniteRange           = 0,	//!< lower and upper limits finite
		kNegativeInfiniteRange = 1,	//!< lower limit infinite, upper limit finite
		kPositiveInfiniteRange = 2,	//!< lower limit finite, upper limit infinite
		kInfiniteRange         = 3,	//!< lower and upper limits infinite
		kEmptyRange            = 4	//!< lower limit == upper limit
	};

	/** @} */
	
	/** \name Constructors and destructors */
	/** @{ */
	
	/**
	 * Constructor */
	BCPrior()
	{ }

	/**
	 * Copy constructor */
	BCPrior(const BCPrior & other)
	{ }

	/**
	 * Destructor */
	virtual ~BCPrior()
	{ }

	/** @} */

	/** \name Getters */
	/** @{ **/

	/** Get as TF1 for drawing purposes.
	 * @param xmin lower limit of range for TF1
	 * @param xmax upper limit of range for TF1
	 * @param normalize whether to normalize TF1 over range*/
	virtual TF1 * GetAsTF1(double xmin=-std::numeric_limits<double>::infinity(), double xmax=std::numeric_limits<double>::infinity(), bool normalize=true) const
	{ return NULL; }

	/**
	 * Get prior
	 * @param x value to evaluate prior at
	 * @return prior */
	virtual double GetPrior(double x) const
	{ return exp(GetLogPrior(x)); }
	
	/**
	 * Get log of prior
	 * @param x value to evaluate log of prior at
	 * @return log of prior */
	virtual double GetLogPrior(double x) const = 0;
	
	/**
	 * Get raw moment of prior distrubion. If limits are infinite, use exact value from prior type.
	 * @param n moment number
	 * @param xmin lower limit of range to evaluate over
	 * @param xmax upper limit of range to evaluate over
	 * @return raw moment of prior distribution */
	virtual double GetRawMoment(unsigned n, double xmin=-std::numeric_limits<double>::infinity(), double xmax=std::numeric_limits<double>::infinity()) const = 0;

	/**
	 * Get integral of prior.
	 * @param xmin lower limit of range to evaluate over
	 * @param xmax upper limit of range to evaluate over
	 * @return integral of prior */
	virtual double GetIntegral(double xmin=-std::numeric_limits<double>::infinity(), double xmax=std::numeric_limits<double>::infinity()) const = 0;

	/**
	 * Get central moment of prior distrubion. If limits are infinite, use exact value from prior type.
	 * @param n moment number
	 * @param xmin lower limit of range to evaluate over
	 * @param xmax upper limit of range to evaluate over
	 * @return central of prior distribution */
	virtual double GetCentralMoment(unsigned n, double xmin=-std::numeric_limits<double>::infinity(), double xmax=std::numeric_limits<double>::infinity()) const;
	
	/**
	 * Get standardised moment of prior distrubion. If limits are infinite, use exact value from prior type.
	 * @param n moment number
	 * @param xmin lower limit of range to evaluate over
	 * @param xmax upper limit of range to evaluate over
	 * @return standardised moment of prior distribution */
	virtual double GetStandardisedMoment(unsigned n, double xmin=-std::numeric_limits<double>::infinity(), double xmax=std::numeric_limits<double>::infinity()) const;

	/**
	 * Get mean of prior. If limits are infinite, use exact value from prior type
	 * @param xmin lower limit of range to evaluate over
	 * @param xmax upper limit of range to evaluate over
	 * @return mean of prior distribution */
	virtual double GetMean(double xmin=-std::numeric_limits<double>::infinity(), double xmax=std::numeric_limits<double>::infinity()) const
	{ return GetRawMoment(1,xmin,xmax); }

	/**
	 * Get variance of prior. If limits are infinite, use exact value from prior type
	 * @param xmin lower limit of range to evaluate over
	 * @param xmax upper limit of range to evaluate over
	 * @return variance of prior distribution */
	virtual double GetVariance(double xmin=-std::numeric_limits<double>::infinity(), double xmax=std::numeric_limits<double>::infinity()) const
	{ return GetCentralMoment(2,xmin,xmax); }

	/**
	 * Get standard deviation of prior. If limits are infinite, use exact value from prior type
	 * @param xmin lower limit of range to evaluate over
	 * @param xmax upper limit of range to evaluate over
	 * @return standard deviation of prior distribution */
	virtual double GetStandardDeviation(double xmin=-std::numeric_limits<double>::infinity(), double xmax=std::numeric_limits<double>::infinity()) const
	{ return sqrt(fabs(GetVariance())); }

	/**
	 * Get skewness of prior. If limits are infinite, use exact value from prior type
	 * @param xmin lower limit of range to evaluate over
	 * @param xmax upper limit of range to evaluate over
	 * @return skewness of prior distribution */
	virtual double GetSkewness(double xmin=-std::numeric_limits<double>::infinity(), double xmax=std::numeric_limits<double>::infinity()) const
	{ return GetStandardisedMoment(3,xmin,xmax); }

	/**
	 * Get kurtosis of prior. If limits are infinite, use exact value from prior type
	 * @param xmin lower limit of range to evaluate over
	 * @param xmax upper limit of range to evaluate over
	 * @return kurtosis of prior distribution */
	virtual double GetKurtosis(double xmin=-std::numeric_limits<double>::infinity(), double xmax=std::numeric_limits<double>::infinity()) const
	{ return GetStandardisedMoment(4,xmin,xmax); }

	/** @} **/

	/** \name Miscellaneous functions. */
	/** @{ **/

	/** Check type of limits
	 * @param xmin lower limit of range
	 * @param xmax upper limit of range
	 * @return type of range limits */
	BCPrior::BCPriorRange CheckLimits(double xmin, double xmax) const;

	/** @} **/
};

/** \class BCSplitGaussianPrior **/
class BCSplitGaussianPrior : public BCPrior {
public:
	BCSplitGaussianPrior(double mean, double sigma_low, double sigma_high) : BCPrior(BCPrior::kPriorSplitGaussian), fMean(mean), fSigmaLow(sigma_low), fSigmaHigh(sigma_high) {}
	virtual ~BCSplitGaussianPrior() {}
	virtual double GetLogPrior(double x) const;
	virtual double GetRawMoment(unsigned n, double xmin=-std::numeric_limits<double>::infinity(), double xmax=std::numeric_limits<double>::infinity()) const;
	virtual double GetIntegral(double xmin=-std::numeric_limits<double>::infinity(), double xmax=std::numeric_limits<double>::infinity()) const;

	double GetSigmaLow() { return fSigmaLow; }
	double GetSigmaHigh() { return fSigmaHigh; }

	void SetMean(double mean) { fMean = mean; }
	void SetSigmaLow(double sigma_low)	{ fSigmaLow = sigma_low; }
	void SetSigmaHigh(double sigma_high)	{ fSigmaLow = sigma_high; }
	void SetParameters(double mean, double sigma_low, double sigma_high)
	{ SetMean(mean); SetSigmaLow(sigma_low); SetSigmaHigh(sigma_high); }

protected:
	double fMean;
	double fSigmaLow;
	double fSigmaHigh;
};

/** \class BCCauchyPrior **/
class BCCauchyPrior : public BCPrior {
public:
	BCCauchyPrior(double mean, double scale) : BCPrior(BCPrior::kPriorCauchy), fMean(mean), fScale(scale) {}
	virtual ~BCCauchyPrior() {}
	virtual double GetLogPrior(double x) const { return -log(1+(x-fMean)*(x-fMean)/fScale/fScale) - log(M_PI*fScale); }
	virtual double GetRawMoment(unsigned n, double xmin=-std::numeric_limits<double>::infinity(), double xmax=std::numeric_limits<double>::infinity()) const;
	virtual double GetIntegral(double xmin=-std::numeric_limits<double>::infinity(), double xmax=std::numeric_limits<double>::infinity()) const;

	void SetMean(double mean) { fMean = mean; }
	void SetScale(double scale)	{ fScale = scale; }
	void SetParameters(double mean, double scale)	{ SetMean(mean); SetScale(scale); }

protected:
	double fMean;
	double fScale;
};

/** \class BCTF1Prior **/
class BCTF1Prior : public BCPrior {
public:
	BCTF1Prior(TF1 * const f) : BCPrior(BCPrior::kPriorTF1), fPriorFunction(new TF1(*f)) {}
	virtual ~BCTF1Prior() {delete fPriorFunction;}

	virtual double GetLogPrior(double x) const
	{return (fPriorFunction) ? log(fPriorFunction->Eval(x)) : -std::numeric_limits<double>::infinity(); }

	virtual double GetRawMoment(unsigned n, double xmin=-std::numeric_limits<double>::infinity(), double xmax=std::numeric_limits<double>::infinity()) const;

	virtual double GetIntegral(double xmin=-std::numeric_limits<double>::infinity(), double xmax=std::numeric_limits<double>::infinity()) const
	{ return (fPriorFunction) ? fPriorFunction->Integral(xmin,xmax) : 0; }

	TF1 * const GetPriorFunction() { return fPriorFunction; }

protected:
	TF1 * fPriorFunction;
};

/** \class BCTH1Prior **/
class BCTH1Prior : public BCPrior {
public:
	BCTH1Prior(TH1 * const h, bool interpolate);
	virtual ~BCTH1Prior() {delete fPriorHistogram;}

	virtual double GetLogPrior(double x) const
	{return (fPriorHistogram) ? log(((fInterpolate) ? fPriorHistogram->Interpolate(x) : fPriorHistogram->GetBinContent(fPriorHistogram->FindFixBin(x)))) : -std::numeric_limits<double>::infinity(); }
	
	virtual double GetRawMoment(unsigned n, double xmin=-std::numeric_limits<double>::infinity(), double xmax=std::numeric_limits<double>::infinity()) const;
	virtual double GetIntegral(double xmin=-std::numeric_limits<double>::infinity(), double xmax=std::numeric_limits<double>::infinity()) const {return 1;}
	virtual double GetCentralMoment(unsigned n, double xmin=-std::numeric_limits<double>::infinity(), double xmax=std::numeric_limits<double>::infinity()) const;
	virtual double GetStandardisedMoment(unsigned n, double xmin=-std::numeric_limits<double>::infinity(), double xmax=std::numeric_limits<double>::infinity()) const;
	virtual double GetMean(double xmin=-std::numeric_limits<double>::infinity(), double xmax=std::numeric_limits<double>::infinity()) const;
	virtual double GetVariance(double xmin=-std::numeric_limits<double>::infinity(), double xmax=std::numeric_limits<double>::infinity()) const;
	virtual double GetStandardDeviation(double xmin=-std::numeric_limits<double>::infinity(), double xmax=std::numeric_limits<double>::infinity()) const;
	virtual double GetSkewness(double xmin=-std::numeric_limits<double>::infinity(), double xmax=std::numeric_limits<double>::infinity()) const;
	virtual double GetKurtosis(double xmin=-std::numeric_limits<double>::infinity(), double xmax=std::numeric_limits<double>::infinity()) const;

	TH1 * const GetPriorHistogram() { return fPriorHistogram; }
	void SetInterpolate(bool interpolate) {fInterpolate=interpolate;}
	bool GetInterpolate() {return fInterpolate;}

protected:
	TH1 * fPriorHistogram;
	bool fInterpolate;
};


#endif
