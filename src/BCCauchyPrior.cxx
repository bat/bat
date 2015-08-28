/*
 * Copyright (C) 2007-2015, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

#include "BCCauchyPrior.h"
#include "BCAux.h"

#include <cmath>

#include <TF1.h>


// ---------------------------------------------------------
BCCauchyPrior::BCCauchyPrior(double mean, double scale)
    : BCPrior()
    , fMean(mean)
    , fScale(scale)
{
}

// ---------------------------------------------------------
double BCCauchyPrior::GetRawMoment(unsigned n, double xmin, double xmax)
{
    if (n == 0)
        return BCPrior::GetRawMoment(n, xmin, xmax);

    BCAux::BCRange r = BCAux::RangeType(xmin, xmax);

    if (r == BCAux::kReverseRange)
        return GetRawMoment(n, xmax, xmin);

    if (r == BCAux::kEmptyRange)
        return (n == 1) ? xmin : 0;

    if (r == BCAux::kInfiniteRange)
        return (n == 1) ? fMean : std::numeric_limits<double>::quiet_NaN();

    if (r == BCAux::kNegativeInfiniteRange)
        return (n == 1) ? -std::numeric_limits<double>::infinity() : std::numeric_limits<double>::quiet_NaN();

    if (r == BCAux::kPositiveInfiniteRange)
        return (n == 1) ? +std::numeric_limits<double>::infinity() : std::numeric_limits<double>::quiet_NaN();

    // else finite

    double L = (xmin - fMean) / fScale;
    double H = (xmax - fMean) / fScale;

    if (n == 1)
        return fMean + fScale / 2 * log((1 + H * H) / (1 + L * L)) / (atan(H) - atan(L));

    if (n == 2)
        return fMean * fMean - fScale * fScale + fScale * (log((1 + H * H) / (1 + L * L)) + xmax - xmin) / (atan(H) - atan(L));

    return BCPrior::GetRawMoment(n, xmin, xmax);
}

// ---------------------------------------------------------
double BCCauchyPrior::GetIntegral(double xmin, double xmax)
{
    switch (BCAux::RangeType(xmin, xmax)) {

        case BCAux::kFiniteRange:
            return (atan((xmax - fMean) / fScale) - atan((xmin - fMean) / fScale)) / M_PI;

        case BCAux::kNegativeInfiniteRange:
            return 0.5 + atan((xmax - fMean) / fScale) / M_PI;

        case BCAux::kPositiveInfiniteRange:
            return 0.5 - atan((xmin - fMean) / fScale) / M_PI;

        case BCAux::kInfiniteRange:
            return 1;

        case BCAux::kEmptyRange:
            return 0;

        default:
            return std::numeric_limits<double>::infinity();
    }
}


