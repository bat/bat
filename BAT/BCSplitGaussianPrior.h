#ifndef __BCSPLITGAUSSIANPRIOR__H
#define __BCSPLITGAUSSIANPRIOR__H

/*!
 * \class BCSplitGaussianPrior
 * \brief A class to represent a split-Gaussian prior of a parameter
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

class BCSplitGaussianPrior : public BCPrior
{
public:
    /** \name Constructor & Destructor */
    /** @{ **/

    /** Constructor */
    BCSplitGaussianPrior(double mode, double sigma_below, double sigma_above);

    /** Copy constructor */
    BCSplitGaussianPrior(const BCSplitGaussianPrior& other);

    /** Destructor */
    virtual ~BCSplitGaussianPrior() {};

    /** @} **/

    /** \name Functions overloaded from BCPrior **/
    /** @{ **/

    /** Clone function */
    virtual BCPrior* Clone() const
    { return new BCSplitGaussianPrior(*this); }

    /**
     * @return Whether everything needed for prior is set and prior can be used. */
    virtual bool IsValid() const
    {
        return std::isfinite(fMode)
               and std::isfinite(fSigmaBelow) and fSigmaBelow > 0
               and std::isfinite(fSigmaAbove) and fSigmaAbove > 0;
    }

    /**
     * Get log of prior
     * @param x value to evaluate log of prior at
     * @return log of prior */
    virtual double GetLogPrior(double x);

    /**
     * Return mode of prior (in range).
     * @param xmin lower limit of range to evaluate over
     * @param xmax upper limit of range to evaluate over
     * @return mode of prior in range. */
    virtual double GetMode(double xmin = -std::numeric_limits<double>::infinity(), double xmax = std::numeric_limits<double>::infinity())
    { if (fMode < xmin) return xmin; if (fMode > xmax) return xmax; return fMode; }

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
    { SetMode(mean); }

    void SetMode(double mode)
    { fMode = mode; }

    void SetSigmaBelow(double sigma)
    { fSigmaBelow = sigma; }

    void SetSigmaAbove(double sigma)
    { fSigmaAbove = sigma; }

    void SetParameters(double mode, double sigma_below, double sigma_above)
    { SetMode(mode); SetSigmaBelow(sigma_below); SetSigmaAbove(sigma_above);}

    /** @} **/

protected:
    double fMode;							 ///< mode of split gaussian
    double fSigmaBelow;				 ///< std dev of split gaussian below mode
    double fSigmaAbove;				 ///< std dev of split gaussian above mode
};

#endif
