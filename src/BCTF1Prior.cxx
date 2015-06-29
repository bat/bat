/*
 * Copyright (C) 2007-2015, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

#include "BCTF1Prior.h"

#include <cmath>

#include <TF1.h>

// ---------------------------------------------------------
BCTF1Prior::BCTF1Prior(TF1& f)
    : BCPrior()
      // fPriorFunction(f)
{
    fPriorFunction = f;
}

// ---------------------------------------------------------
BCTF1Prior::BCTF1Prior(const char* formula, double xmin, double xmax)
    : BCPrior()
{
    fPriorFunction = TF1("f1_prior", formula, xmin, xmax);
}

// ---------------------------------------------------------
BCTF1Prior::BCTF1Prior(const BCTF1Prior& other)
    : BCPrior(other)
{
    fPriorFunction = other.fPriorFunction;
}

// ---------------------------------------------------------
BCTF1Prior::~BCTF1Prior()
{
}






