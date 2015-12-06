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
BCTF1LogPrior::BCTF1LogPrior(TF1& f)
    : BCPrior(),
      fLogPriorFunction(f)
{
}

// ---------------------------------------------------------
BCTF1LogPrior::BCTF1LogPrior(const std::string& formula, double xmin, double xmax)
    : BCPrior(),
      fLogPriorFunction("f1_logprior", formula.c_str(), xmin, xmax)
{
}

// ---------------------------------------------------------
BCTF1LogPrior::~BCTF1LogPrior()
{
}

// ---------------------------------------------------------
void BCTF1LogPrior::SetFunctionRange(double xmin, double xmax)
{
    BCPrior::SetFunctionRange(xmin, xmax);
    fLogPriorFunction.SetRange(xmin, xmax);
}

// ---------------------------------------------------------
double BCTF1LogPrior::GetMode(double xmin, double xmax)
{
    BCAux::MakeFinite(xmin, xmax);
    return fLogPriorFunction.GetMaximumX(xmin, xmax);
}
