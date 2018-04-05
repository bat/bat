/*
 * Copyright (C) 2007-2018, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

#include "BCPositiveDefinitePrior.h"

// ---------------------------------------------------------
BCPositiveDefinitePrior::BCPositiveDefinitePrior(BCPrior* prior)
    : BCPrior(),
      fPrior(prior)
{
}

// ---------------------------------------------------------
BCPositiveDefinitePrior::BCPositiveDefinitePrior(const BCPositiveDefinitePrior& other)
    : BCPrior(other),
      fPrior(other.fPrior->Clone())
{
}

// ---------------------------------------------------------
BCPositiveDefinitePrior::~BCPositiveDefinitePrior()
{
    delete fPrior;
}

// ---------------------------------------------------------
BCPositiveDefinitePrior& BCPositiveDefinitePrior::operator=(BCPositiveDefinitePrior other)
{
    swap(*this, other);
    return *this;
}

// ---------------------------------------------------------
void swap(BCPositiveDefinitePrior& A, BCPositiveDefinitePrior& B)
{
    swap(static_cast<BCPrior&>(A), static_cast<BCPrior&>(B));
    std::swap(A.fPrior, B.fPrior);
}

// ---------------------------------------------------------
double BCPositiveDefinitePrior::GetRandomValue(double xmin, double xmax, TRandom* const R)
{
    xmin = std::max<double>(xmin, 0);
    xmax = std::max<double>(xmax, 0);
    if (xmin == xmax)
        return std::numeric_limits<double>::quiet_NaN();
    return fPrior->GetRandomValue(xmin, xmax, R);
}
