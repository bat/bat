/*
 * Copyright (C) 2007-2018, the BAT core developer team
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
{
    fPriorFunction = f;
}

// ---------------------------------------------------------
BCTF1Prior::BCTF1Prior(const std::string& formula, double xmin, double xmax)
    : BCPrior()
{
    fPriorFunction = TF1("f1_prior", formula.c_str(), xmin, xmax);
}
