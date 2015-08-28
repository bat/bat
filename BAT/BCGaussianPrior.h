#ifndef __BCGAUSSIANPRIOR__H
#define __BCGAUSSIANPRIOR__H

/*!
 * \class BCGaussianPrior
 * \brief A class to represent a gaussian prior of a parameter
 * \author Daniel Greenwald
 * \version 1.0
 * \date 01.2015
 * \ingroup Priors
 */


/*
 * Copyright (C) 2007-2015, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

// ---------------------------------------------------------

#include "BCPrior.h"

#include <limits>
#include <cmath>

class TF1;

// ---------------------------------------------------------

class BCGaussianPrior : public BCPrior
{
public:
    /** \name Constructor & Destructor */
    /** @{ **/

    /** Constructor */
    BCGaussianPrior(double mean, double sigma);

    /** Destructor */
    virtual ~BCGaussianPrior() {};

    /** @} **/

    /** \name Functions overloaded from BCPrior **/
    /** @{ **/

    /** Clone function */
    virtual BCPrior* Clone() const
    { return new BCGaussianPrior(*this); }

    /**
     * @return Whether everything needed for prior is set and prior can be used. */
    virtual bool IsValid() const
    { return std::isfinite(fMean) and std::isfinite(fSigma) and fSigma > 0;}

    /**
     * Get log of prior
     * @param x value to evaluate log of prior at
     * @return log of prior */
    virtual double GetLogPrior(double x)
    { return -0.5 * (x - fMean) * (x - fMean) / fSigma / fSigma - log(fSigma) - 0.5 * log(2 * M_PI); }

    /**
     * Return mode of prior (in range).
     * @param xmin lower limit of range to evaluate over
     * @param xmax upper limit of range to evaluate over
     * @return mode of prior in range. */
    virtual double GetMode(double xmin = -std::numeric_limits<double>::infinity(), double xmax = std::numeric_limits<double>::infinity())
    { if (fMean < xmin) return xmin; if (fMean > xmax) return xmax; return fMean; }

    /**
     * Get raw moment of prior distrubion. If limits are infinite, use exact value from prior type.
     * @param n moment number
     * @param xmin lower limit of range to evaluate over
     * @param xmax upper limit of range to evaluate over
     * @return raw moment of prior distribution */
    virtual double GetRawMoment(unsigned n, double xmin = -std::numeric_limits<double>::infinity(), double xmax = std::numeric_limits<double>::infinity());

    /**
     * Get integral of prior.
     * @param xmin lower limit of range to evaluate over
     * @param xmax upper limit of range to evaluate over
     * @return integral of prior */
    virtual double GetIntegral(double xmin = -std::numeric_limits<double>::infinity(), double xmax = std::numeric_limits<double>::infinity());

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

#endif
