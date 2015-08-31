#ifndef __BCFORMULALOGPRIOR__H
#define __BCFORMULALOGPRIOR__H

/*!
 * \class BCPriorLogTF1
 * \brief A class to represent the log of a prior of a parameter by a formula through a TF1
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

#include <TF1.h>

// ---------------------------------------------------------

class BCPriorLogTF1 : public BCPrior
{
public:
    /** \name Constructor & Destructor */
    /** @{ **/

    /** Constructor taking a TF1**/
    BCPriorLogTF1(TF1& f);

    /** Constructor taking a formula**/
    BCPriorLogTF1(const char* formula, double xmin, double xmax); //double xmin=-std::numeric_limits<double>::infinity(),xmax=std::numeric_limits<double>::infinity());

    /** Destructor */
    virtual ~BCPriorLogTF1();

    /** @} **/


    /** \name Functions overloaded from BCPrior **/
    /** @{ **/

    /** Clone function */
    virtual BCPrior* Clone() const
    { return new BCPriorLogTF1(*this); }

    /**
     * @return Whether everything needed for prior is set and prior can be used. */
    virtual bool IsValid() const
    { return true; }

    /** Set ROOT function range. */
    virtual void SetFunctionRange(double xmin, double xmax);

    /**
     * Get log of prior
     * @param x value to evaluate log of prior at
     * @return log of prior */
    virtual double GetLogPrior(double x)
    { return fLogPriorFunction.Eval(x); }

    /**
     * Return mode of prior (in range).
     * @param xmin lower limit of range to evaluate over
     * @param xmax upper limit of range to evaluate over
     * @return mode of prior in range. */
    virtual double GetMode(double xmin = -std::numeric_limits<double>::infinity(), double xmax = std::numeric_limits<double>::infinity());

    /** @} **/

    /** \name Getters */
    /** @{ **/

    TF1& GetLogFunction()
    { return fLogPriorFunction; }

    const TF1& GetLogFunction() const
    { return fLogPriorFunction; }

    /** @} **/

protected:
    TF1 fLogPriorFunction;			//< TF1 holding log(prior)
};

#endif
