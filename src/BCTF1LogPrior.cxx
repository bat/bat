/*
 * Copyright (C) 2007-2015, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

#include "BCTF1LogPrior.h"
#include "BCAux.h"

#include <cmath>

// ---------------------------------------------------------
BCTF1LogPrior::BCTF1LogPrior(const TF1* const f)
    : BCPrior()
    , fLogPriorFunction(NULL)
{
    if (f)
        fLogPriorFunction = new TF1(*f);
}


// ---------------------------------------------------------
BCTF1LogPrior::BCTF1LogPrior(const char* formula, double xmin, double xmax)
    : BCPrior()
    , fLogPriorFunction(new TF1("f1_logprior", formula, xmin, xmax))
{
}

// ---------------------------------------------------------
BCTF1LogPrior::BCTF1LogPrior(const BCTF1LogPrior& other)
    : BCPrior(other)
{
    if (other.fLogPriorFunction)
        fLogPriorFunction = new TF1(*fLogPriorFunction);
}

// ---------------------------------------------------------
BCTF1LogPrior::~BCTF1LogPrior()
{
    delete fLogPriorFunction;
}

// ---------------------------------------------------------
void BCTF1LogPrior::SetFunctionRange(double xmin, double xmax)
{
    BCPrior::SetFunctionRange(xmin, xmax);
    if (fLogPriorFunction)
        fLogPriorFunction -> SetRange(xmin, xmax);
}

// ---------------------------------------------------------
double BCTF1LogPrior::GetMode(double xmin, double xmax) const
{
    if (!fLogPriorFunction)
        return std::numeric_limits<double>::quiet_NaN();
    BCAux::MakeFinite(xmin, xmax);
    return fLogPriorFunction->GetMaximumX(xmin, xmax);
}
