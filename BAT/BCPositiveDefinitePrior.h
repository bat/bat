#ifndef __BCPOSITIVEDEFINITEPRIOR__H
#define __BCPOSITIVEDEFINITEPRIOR__H

/*!
 * \class BCPositiveDefinitePrior
 * \brief A class to wrap around a BCPrior to make it positive definite
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

#include <algorithm>
#include <limits>

#include "BCPrior.h"

class TRandom;

// ---------------------------------------------------------

class BCPositiveDefinitePrior : public BCPrior
{
public:

    /** \name Constructors and destructors */
    /** @{ */

    /**
     * Constructor */
    BCPositiveDefinitePrior(BCPrior* prior);

    /**
     * Copy constructor */
    BCPositiveDefinitePrior(const BCPositiveDefinitePrior& other);

    /**
     * Destructor */
    virtual ~BCPositiveDefinitePrior();

    /** @} */

    /** \name operator and swap */
    /** @{ **/

    /** assignment operator */
    BCPositiveDefinitePrior& operator=(BCPositiveDefinitePrior other);

    /** swap */
    friend void swap(BCPositiveDefinitePrior& A, BCPositiveDefinitePrior& B);

    /** @} **/

    /** \name Access to wrapped prior. */
    /** @{ */

    BCPrior* Prior()
    { return fPrior; }

    const BCPrior* Prior() const
    { return fPrior; }

    /** @} */

    /** \name Functions overloaded from BCPrior */
    /** @{ */

    /**
     * Get log of prior
     * @param x value to evaluate log of prior at
     * @return log of prior */
    virtual double GetLogPrior(double x)
    { return (x >= 0 and fPrior) ? fPrior->GetLogPrior(x) : -std::numeric_limits<double>::infinity(); }

    /**
     * Get prior
     * @param x value to evaluate prior at
     * @param normalize whether to normalize with stored integral
     * @return prior */
    virtual double GetPrior(double x, bool normalize = false)
    { return (x >= 0 and fPrior) ? fPrior->GetPrior(x, normalize) : 0; }

    /**
     * Clone function */
    virtual BCPrior* Clone() const
    { return new BCPositiveDefinitePrior(*this); }

    /**
     * @return Whether everything needed for prior is set and prior can be used. */
    virtual bool IsValid() const
    { return fPrior != NULL and fPrior->IsValid(); }

    /** Set ROOT function range. */
    virtual void SetFunctionRange(double xmin, double xmax)
    { BCPrior::SetFunctionRange(xmin, xmax); if (fPrior) fPrior->SetFunctionRange(xmin, xmax); }

    /**
     * Return mode of prior (in range).
     * @param xmin lower limit of range to evaluate over
     * @param xmax upper limit of range to evaluate over
     * @return mode of prior in range. */
    virtual double GetMode(double xmin = -std::numeric_limits<double>::infinity(), double xmax = std::numeric_limits<double>::infinity())
    { return (fPrior) ? fPrior->GetMode(std::max<double>(xmin, 0), std::max<double>(xmax, 0)) : std::numeric_limits<double>::quiet_NaN(); }

    /**
     * Get raw moment of prior distrubion. If limits are infinite, use exact value from prior type.
     * @param n moment number
     * @param xmin lower limit of range to evaluate over
     * @param xmax upper limit of range to evaluate over
     * @return raw moment of prior distribution */
    virtual double GetRawMoment(unsigned n, double xmin = -std::numeric_limits<double>::infinity(), double xmax = std::numeric_limits<double>::infinity())
    { return (fPrior) ? fPrior->GetRawMoment(n, std::max<double>(xmin, 0), std::max<double>(xmax, 0)) : std::numeric_limits<double>::quiet_NaN(); }

    /**
     * Get integral of prior.
     * @param xmin lower limit of range to evaluate over
     * @param xmax upper limit of range to evaluate over
     * @return integral of prior */
    virtual double GetIntegral(double xmin = -std::numeric_limits<double>::infinity(), double xmax = std::numeric_limits<double>::infinity())
    { return (fPrior) ? fPrior->GetIntegral(std::max<double>(xmin, 0), std::max<double>(xmax, 0)) : std::numeric_limits<double>::quiet_NaN(); }

    /**
     * Get central moment of prior distrubion. If limits are infinite, use exact value from prior type.
     * @param n moment number
     * @param xmin lower limit of range to evaluate over
     * @param xmax upper limit of range to evaluate over
     * @return central of prior distribution */
    virtual double GetCentralMoment(unsigned n, double xmin = -std::numeric_limits<double>::infinity(), double xmax = std::numeric_limits<double>::infinity())
    { return (fPrior) ? fPrior->GetCentralMoment(n, std::max<double>(xmin, 0), std::max<double>(xmax, 0)) : std::numeric_limits<double>::quiet_NaN(); }

    /**
     * Get standardised moment of prior distrubion. If limits are infinite, use exact value from prior type.
     * @param n moment number
     * @param xmin lower limit of range to evaluate over
     * @param xmax upper limit of range to evaluate over
     * @return standardised moment of prior distribution */
    virtual double GetStandardizedMoment(unsigned n, double xmin = -std::numeric_limits<double>::infinity(), double xmax = std::numeric_limits<double>::infinity())
    { return (fPrior) ? fPrior->GetCentralMoment(n, std::max<double>(xmin, 0), std::max<double>(xmax, 0)) : std::numeric_limits<double>::quiet_NaN(); }

    /**
     * Get mean of prior. If limits are infinite, use exact value from prior type
     * @param xmin lower limit of range to evaluate over
     * @param xmax upper limit of range to evaluate over
     * @return mean of prior distribution */
    virtual double GetMean(double xmin = -std::numeric_limits<double>::infinity(), double xmax = std::numeric_limits<double>::infinity())
    { return (fPrior) ? fPrior->GetMean(std::max<double>(xmin, 0), std::max<double>(xmax, 0)) : std::numeric_limits<double>::quiet_NaN(); }

    /**
     * Get variance of prior. If limits are infinite, use exact value from prior type
     * @param xmin lower limit of range to evaluate over
     * @param xmax upper limit of range to evaluate over
     * @return variance of prior distribution */
    virtual double GetVariance(double xmin = -std::numeric_limits<double>::infinity(), double xmax = std::numeric_limits<double>::infinity())
    { return (fPrior) ? fPrior->GetVariance(std::max<double>(xmin, 0), std::max<double>(xmax, 0)) : std::numeric_limits<double>::quiet_NaN(); }

    /**
     * Get standard deviation of prior. If limits are infinite, use exact value from prior type
     * @param xmin lower limit of range to evaluate over
     * @param xmax upper limit of range to evaluate over
     * @return standard deviation of prior distribution */
    virtual double GetStandardDeviation(double xmin = -std::numeric_limits<double>::infinity(), double xmax = std::numeric_limits<double>::infinity())
    { return (fPrior) ? fPrior->GetStandardDeviation(std::max<double>(xmin, 0), std::max<double>(xmax, 0)) : std::numeric_limits<double>::quiet_NaN(); }

    /**
     * Get skewness of prior. If limits are infinite, use exact value from prior type
     * @param xmin lower limit of range to evaluate over
     * @param xmax upper limit of range to evaluate over
     * @return skewness of prior distribution */
    virtual double GetSkewness(double xmin = -std::numeric_limits<double>::infinity(), double xmax = std::numeric_limits<double>::infinity())
    { return (fPrior) ? fPrior->GetSkewness(std::max<double>(xmin, 0), std::max<double>(xmax, 0)) : std::numeric_limits<double>::quiet_NaN(); }

    /**
     * Get kurtosis of prior. If limits are infinite, use exact value from prior type
     * @param xmin lower limit of range to evaluate over
     * @param xmax upper limit of range to evaluate over
     * @return kurtosis of prior distribution */
    virtual double GetKurtosis(double xmin = -std::numeric_limits<double>::infinity(), double xmax = std::numeric_limits<double>::infinity())
    { return (fPrior) ? fPrior->GetKurtosis(std::max<double>(xmin, 0), std::max<double>(xmax, 0)) : std::numeric_limits<double>::quiet_NaN(); }

    /**
     * @return a random value distributed according to the prior.
     * @param xmin lower limit of range to generate value in
     * @param xmax upper limit of range to generate value in
     * @param R Pointer to the random generator to be used, if needed.
     * @return random value. */
    virtual double GetRandomValue(double xmin, double xmax, TRandom* const R = NULL);

    /** @} */

protected:

    BCPrior* fPrior;

};

#endif
