#ifndef __BCFORMULALOGPRIOR__H
#define __BCFORMULALOGPRIOR__H

/*!
 * \class BCTF1LogPrior
 * \brief A class to represent the log of a prior of a parameter by a formula through a TF1
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

#include "BCPrior.h"

#include <limits>

#include <TF1.h>

// ---------------------------------------------------------

class BCTF1LogPrior : public BCPrior
{
public:
    /** \name Constructor & Destructor */
    /** @{ **/

    /** Constructor taking a TF1**/
    BCTF1LogPrior(const TF1* const f);

    /** Constructor taking a formula**/
    BCTF1LogPrior(const char* formula, double xmin, double xmax); //double xmin=-std::numeric_limits<double>::infinity(),xmax=std::numeric_limits<double>::infinity());

    /** Copy constructor */
    BCTF1LogPrior(const BCTF1LogPrior& other);

    /** Destrcutor */
    virtual ~BCTF1LogPrior();

    /** @} **/

    /** \name Functions overloaded from BCPrior **/
    /** @{ **/

    /** Clone function */
    virtual BCPrior* Clone() const
    { return new BCTF1LogPrior(*this); }

    /**
     * @return Whether everything needed for prior is set and prior can be used. */
    virtual bool IsValid() const
    { return fLogPriorFunction != NULL; }

    /** Set ROOT function range. */
    virtual void SetFunctionRange(double xmin, double xmax);

    /**
     * Get log of prior
     * @param x value to evaluate log of prior at
     * @return log of prior */
    virtual double GetLogPrior(double x) const
    { return (fLogPriorFunction) ? fLogPriorFunction->Eval(x) : -std::numeric_limits<double>::infinity(); }

    /**
     * Return mode of prior (in range).
     * @param xmin lower limit of range to evaluate over
     * @param xmax upper limit of range to evaluate over
     * @return mode of prior in range. */
    virtual double GetMode(double xmin = -std::numeric_limits<double>::infinity(), double xmax = std::numeric_limits<double>::infinity()) const;

    /** @} **/

    /** \name Getters */
    /** @{ **/

    TF1* GetLogFunction()
    { return fLogPriorFunction; }

    /** @} **/

protected:
    TF1* fLogPriorFunction;			//< TF1 holding log(prior)
};

#endif
