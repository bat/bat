#ifndef __BCCAUCHYPRIOR__H
#define __BCCAUCHYPRIOR__H

/*!
 * \class BCCauchyPrior
 * \brief A class to represent a Cauchy prior of a parameter
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

class TF1;

// ---------------------------------------------------------

class BCCauchyPrior : public BCPrior
{
public:
    /** \name Constructor & Destructor */
    /** @{ **/

    /** Constructor */
    BCCauchyPrior(double mean, double scale);

    /** Destrcutor */
    virtual ~BCCauchyPrior() {};

    /** @} **/

    /** \name Functions overloaded from BCPrior **/
    /** @{ **/

    /** Clone function */
    virtual BCPrior* Clone() const
    { return new BCCauchyPrior(*this); }

    /**
     * @return Whether everything needed for prior is set and prior can be used. */
    virtual bool IsValid() const
    { return std::isfinite(fMean) and std::isfinite(fScale) and fScale > 0;}

    /**
     * Get log of prior
     * @param x value to evaluate log of prior at
     * @return log of prior */
    virtual double GetLogPrior(double x)
    { return log(fScale) - log(fScale * fScale + (x - fMean) * (x - fMean)) - log(M_PI); }

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

    virtual double GetScale() const
    { return fScale; }

    /** @} **/

    /** \name Setters **/
    /** @{ **/

    virtual void SetMean(double mean)
    { fMean = mean; }

    virtual void SetScale(double scale)
    { fScale = scale; }

    virtual void SetParameters(double mean, double scale)
    { SetMean(mean); SetScale(scale); }

    /** @} **/

protected:
    double fMean;							 ///< mean of Cauchy distribution
    double fScale;						 ///< scale of Cauchy distribution
};

#endif
