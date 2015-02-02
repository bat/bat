#ifndef __BCGAUSSIANPRIOR__H
#define __BCGAUSSIANPRIOR__H

/*!
 * \class BCGaussianPrior
 * \brief A class to represent a gaussian prior of a parameter
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

class BCGaussianPrior : public BCPrior {
public:
	/** \name Constructor & Destructor */
	/** @{ **/

	/** Constructor */
	BCGaussianPrior(double mean, double sigma);

	/** Copy constructor */
	BCGaussianPrior(const BCGaussianPrior & other);

	/** Destrcutor */
	virtual ~BCGaussianPrior() {}

	/** @} **/

	/** \name Functions overloaded from BCPrior **/
	/** @{ **/

	/** Get as TF1 for drawing purposes.
	 * Parameter zero is mean, parameter one is sigma,
	 * last parameter is integral over range
	 * @param xmin lower limit of range for TF1
	 * @param xmax upper limit of range for TF1
	 * @param normalize whether to normalize TF1 over range*/
	virtual TF1 * GetAsTF1(double xmin=-std::numeric_limits<double>::infinity(), double xmax=std::numeric_limits<double>::infinity(), bool normalize=true) const;

	/**
	 * Get log of prior
	 * @param x value to evaluate log of prior at
	 * @return log of prior */
	virtual double GetLogPrior(double x) const { return -(x-fMean)*(x-fMean)/fSigma/fSigma/2. - log(fSigma) - log(2*M_PI); }

	/**
	 * Get raw moment of prior distrubion. If limits are infinite, use exact value from prior type.
	 * @param n moment number
	 * @param xmin lower limit of range to evaluate over
	 * @param xmax upper limit of range to evaluate over
	 * @return raw moment of prior distribution */
	virtual double GetRawMoment(unsigned n, double xmin=-std::numeric_limits<double>::infinity(), double xmax=std::numeric_limits<double>::infinity()) const;

	/**
	 * Get integral of prior.
	 * @param xmin lower limit of range to evaluate over
	 * @param xmax upper limit of range to evaluate over
	 * @return integral of prior */
	virtual double GetIntegral(double xmin=-std::numeric_limits<double>::infinity(), double xmax=std::numeric_limits<double>::infinity()) const;

	/** @} **/

	/** \name Setters **/
	/** @{ **/

	void SetMean(double mean)
	{ fMean = mean; }

	void SetSigma(double sigma)
	{ fSigma = sigma; }

	void SetParameters(double mean, double sigma)
	{ SetMean(mean); SetSigma(sigma); }

	/** @} **/

protected:
	double fMean;									///< mean of Gaussian
	double fSigma;								///< std dev of Gaussian
};
