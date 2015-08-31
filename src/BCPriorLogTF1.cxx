/*
 * Copyright (C) 2007-2015, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

#include "BCPriorLogTF1.h"

#include "BCAux.h"

#include <cmath>

// ---------------------------------------------------------
BCPriorLogTF1::BCPriorLogTF1(TF1& f)
    : BCPrior(),
      fLogPriorFunction(f)
{
}

// ---------------------------------------------------------
BCPriorLogTF1::BCPriorLogTF1(const char* formula, double xmin, double xmax)
    : BCPrior(),
      fLogPriorFunction("f1_logprior", formula, xmin, xmax)
{
}

// ---------------------------------------------------------
BCPriorLogTF1::~BCPriorLogTF1()
{
}

// ---------------------------------------------------------
void BCPriorLogTF1::SetFunctionRange(double xmin, double xmax)
{
    BCPrior::SetFunctionRange(xmin, xmax);
    fLogPriorFunction.SetRange(xmin, xmax);
}

// ---------------------------------------------------------
double BCPriorLogTF1::GetMode(double xmin, double xmax)
{
    BCAux::MakeFinite(xmin, xmax);
    return fLogPriorFunction.GetMaximumX(xmin, xmax);
}
