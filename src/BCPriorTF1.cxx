/*
 * Copyright (C) 2007-2015, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

#include "BCPriorTF1.h"

#include <cmath>

#include <TF1.h>

// ---------------------------------------------------------
BCPriorTF1::BCPriorTF1(TF1& f)
    : BCPrior()
{
    fPriorFunction = f;
}

// ---------------------------------------------------------
BCPriorTF1::BCPriorTF1(const char* formula, double xmin, double xmax)
    : BCPrior()
{
    fPriorFunction = TF1("f1_prior", formula, xmin, xmax);
}
