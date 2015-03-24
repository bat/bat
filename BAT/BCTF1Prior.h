#ifndef __BCTF1PRIOR__H
#define __BCTF1PRIOR__H

/*!
 * \class BCTF1Prior
 * \brief A class to represent the prior of a parameter by a formula through a TF1
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

class BCTF1Prior : public BCPrior
{
public:
    /** \name Constructor & Destructor */
    /** @{ **/

    /** Constructor taking TF1*/
    BCTF1Prior(const TF1* const f);

    /** Constructor with formula and limits. */
    BCTF1Prior(const char* formula, double xmin, double xmax); //double xmin=-std::numeric_limits<double>::infinity(), double xmax=std::numeric_limits<double>::infinity());

    /** Copy constructor */
    BCTF1Prior(const BCTF1Prior& other);

    /** Destrcutor */
    virtual ~BCTF1Prior();

    /** @} **/

    /** \name Functions overloaded from BCPrior **/
    /** @{ **/

    /** Clone function */
    virtual BCPrior* Clone() const
    { return new BCTF1Prior(*this); }

    /**
     * @return Whether everything needed for prior is set and prior can be used. */
    virtual bool IsValid() const
    { return fPriorFunction != NULL; }

    /**
     * Get prior
     * @param x value to evaluate prior at
     * @return prior */
    virtual double GetPrior(double x) const
    { return (fPriorFunction) ? fPriorFunction->Eval(x) : 0; }

    /**
     * Get prior
     * @param x value to evaluate log(prior) at
     * @return log(prior) */
    virtual double GetLogPrior(double x) const
    { return (fPriorFunction) ? log(fPriorFunction->Eval(x)) : -std::numeric_limits<double>::infinity(); }


    /** @} **/

};

#endif
