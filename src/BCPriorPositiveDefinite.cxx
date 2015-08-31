/*
 * Copyright (C) 2007-2015, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

#include "BCPriorPositiveDefinite.h"

// ---------------------------------------------------------
BCPriorPositiveDefinite::BCPriorPositiveDefinite(BCPrior* prior)
    : BCPrior(),
      fPrior(prior)
{
}

// ---------------------------------------------------------
BCPriorPositiveDefinite::BCPriorPositiveDefinite(const BCPriorPositiveDefinite& other)
    : BCPrior(other),
      fPrior(other.fPrior->Clone())
{
}

// ---------------------------------------------------------
BCPriorPositiveDefinite::~BCPriorPositiveDefinite()
{
    delete fPrior;
}

// ---------------------------------------------------------
BCPriorPositiveDefinite& BCPriorPositiveDefinite::operator=(BCPriorPositiveDefinite other)
{
    swap(*this, other);
    return *this;
}

// ---------------------------------------------------------
void swap(BCPriorPositiveDefinite& A, BCPriorPositiveDefinite& B)
{
    swap(static_cast<BCPrior&>(A), static_cast<BCPrior&>(B));
    std::swap(A.fPrior, B.fPrior);
}

// ---------------------------------------------------------
double BCPriorPositiveDefinite::GetRandomValue(double xmin, double xmax, TRandom* const R)
{
    xmin = std::max<double>(xmin, 0);
    xmax = std::max<double>(xmax, 0);
    if (xmin == xmax)
        return std::numeric_limits<double>::quiet_NaN();
    return fPrior->GetRandomValue(xmin, xmax, R);
}
