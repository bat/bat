#ifndef __BCCONSTANTPRIOR__H
#define __BCCONSTANTPRIOR__H

/*!
 * \class BCConstantPrior
 * \brief A class to represent a constant prior of a parameter
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
class TRandom;

// ---------------------------------------------------------

class BCConstantPrior : public BCPrior
{
public:

    /** \name Constructors and destructor */
    /** @{ */

    /** Constructor for constant unit prior */
    BCConstantPrior();

    /** Constructor for constant 1/range prior */
    BCConstantPrior(double range_width);

    /** Constructor for constant 1/range prior */
    BCConstantPrior(double xmin, double xmax);

    /** Destructor */
    virtual ~BCConstantPrior() {};

    /** @} */

    /** Clone function */
    virtual BCPrior* Clone() const
    { return new BCConstantPrior(*this); }

    /** @return constant log(prior) */
    virtual double GetLogPrior(double /*x*/)
    { return -fLogRangeWidth; }

    /**
     * @return Whether everything needed for prior is set and prior can be used. */
    virtual bool IsValid() const
    { return std::isfinite(fLogRangeWidth); }

    /**
     * Return mode of prior (in range) --- center of interval for constant prior.
     * @param xmin lower limit of range to evaluate over
     * @param xmax upper limit of range to evaluate over
     * @return mode of prior in range. */
    virtual double GetMode(double xmin = -std::numeric_limits<double>::infinity(), double xmax = std::numeric_limits<double>::infinity());

    /** @return raw moments of uniform continuous distribtion. */
    virtual double GetRawMoment(unsigned n, double xmin = -std::numeric_limits<double>::infinity(), double xmax = std::numeric_limits<double>::infinity());

    /** @return integral = 1 */
    virtual double GetIntegral(double xmin = -std::numeric_limits<double>::infinity(), double xmax = std::numeric_limits<double>::infinity())
    { (void)xmin; (void)xmax; return 1; }

    /**
     * @return a random value distributed according to the prior.
     * @param xmin lower limit of range to generate value in
     * @param xmax upper limit of range to generate value in
     * @param R Pointer to the random generator to be used, if needed.
     * @return random value. */
    virtual double GetRandomValue(double xmin, double xmax, TRandom* const R = NULL);

protected:
    ///> log of the width of the parameter range
    double fLogRangeWidth;

};

#endif
