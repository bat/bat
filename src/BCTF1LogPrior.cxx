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
BCTF1LogPrior::BCTF1LogPrior(const char* formula, double xmin, double xmax)
    : BCPrior(),
      fLogPriorFunction("f1_logprior", formula, xmin, xmax)
{
}

// ---------------------------------------------------------
BCTF1LogPrior::BCTF1LogPrior(const BCTF1LogPrior& other)
    : BCPrior(other),
      fLogPriorFunction(other.fLogPriorFunction)
{
}

// ---------------------------------------------------------
BCTF1LogPrior::~BCTF1LogPrior()
{
}

// ---------------------------------------------------------
BCTF1LogPrior& BCTF1LogPrior::operator=(const BCTF1LogPrior& rhs)
{
    BCTF1LogPrior temp(rhs);
    swap(*this, temp);
    return *this;
}

// ---------------------------------------------------------
void swap(BCTF1LogPrior& A, BCTF1LogPrior& B)
{
    swap(static_cast<BCPrior&>(A), static_cast<BCPrior&>(B));
    TF1 temp = A.fLogPriorFunction;
    A.fLogPriorFunction = B.fLogPriorFunction;
    B.fLogPriorFunction = temp;
}

// ---------------------------------------------------------
void BCTF1LogPrior::SetFunctionRange(double xmin, double xmax)
{
    BCPrior::SetFunctionRange(xmin, xmax);
    fLogPriorFunction.SetRange(xmin, xmax);
}

// ---------------------------------------------------------
double BCTF1LogPrior::GetMode(double xmin, double xmax) const
{
    BCAux::MakeFinite(xmin, xmax);
    return fLogPriorFunction.GetMaximumX(xmin, xmax);
}
