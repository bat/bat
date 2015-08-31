#ifndef __BCTF1PRIOR__H
#define __BCTF1PRIOR__H

/*!
 * \class BCPriorTF1
 * \brief A class to represent the prior of a parameter by a formula through a TF1
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

class BCPriorTF1 : public BCPrior
{
public:
    /** \name Constructor & Destructor */
    /** @{ **/

    /** Constructor taking TF1*/
    BCPriorTF1(TF1& f);

    /** Constructor with formula and limits. */
    BCPriorTF1(const char* formula, double xmin, double xmax); //double xmin=-std::numeric_limits<double>::infinity(), double xmax=std::numeric_limits<double>::infinity());

    /** Destrcutor */
    virtual ~BCPriorTF1() {};

    /** @} **/

    /** \name Functions overloaded from BCPrior **/
    /** @{ **/

    /** Clone function */
    virtual BCPrior* Clone() const
    { return new BCPriorTF1(*this); }

    /**
     * @return Whether everything needed for prior is set and prior can be used. */
    virtual bool IsValid() const
    { return true; }

    /**
     * Get prior
     * @param x value to evaluate prior at
     * @param normalize Whether to normalize prior with stored integral
     * @return prior */
    virtual double GetPrior(double x, bool normalize = false)
    { return fPriorFunction.Eval(x) * ((normalize) ? exp(-fLogIntegral) : 1); }

    /**
     * Get prior
     * @param x value to evaluate log(prior) at
     * @return log(prior) */
    virtual double GetLogPrior(double x)
    { return log(fPriorFunction.Eval(x)); }

    /** @} **/

};

#endif
