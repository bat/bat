/*
 * Copyright (C) 2007-2015, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

#include "BCConstantPrior.h"
#include "BCAux.h"

#include <TRandom.h>

// ---------------------------------------------------------
BCConstantPrior::BCConstantPrior()
    : BCPrior()
    , fLogRangeWidth(0)
{}

// ---------------------------------------------------------
BCConstantPrior::BCConstantPrior(double range_width)
    : BCPrior()
    , fLogRangeWidth(0)
{
    if (range_width > 0)
        fLogRangeWidth = log(range_width);
}

// ---------------------------------------------------------
BCConstantPrior::BCConstantPrior(double xmin, double xmax)
    : BCPrior()
    , fLogRangeWidth(0)
{
    if ((xmax - xmin) > 0)
        fLogRangeWidth = log(xmax - xmin);
}

// ---------------------------------------------------------
double BCConstantPrior::GetMode(double xmin, double xmax)
{
    BCAux::BCRange r = BCAux::RangeType(xmin, xmax);

    if (r == BCAux::kReverseRange)
        return GetMode(xmax, xmin);

    if (r == BCAux::kFiniteRange or r == BCAux::kEmptyRange)
        return 0.5 * (xmin + xmax);

    if (r == BCAux::kInfiniteRange)
        return 0;

    if (r == BCAux::kPositiveInfiniteRange)
        return std::numeric_limits<double>::infinity();

    if (r == BCAux::kNegativeInfiniteRange)
        return -std::numeric_limits<double>::infinity();

    return std::numeric_limits<double>::quiet_NaN();
}

// ---------------------------------------------------------
double BCConstantPrior::GetRawMoment(unsigned n, double xmin, double xmax)
{
    if (n == 0)
        return 1;

    BCAux::BCRange r = BCAux::RangeType(xmin, xmax);

    if (r == BCAux::kReverseRange)
        return GetRawMoment(n, xmax, xmin);

    if (r == BCAux::kEmptyRange)
        return (n == 1) ? xmin : 0;

    if (r == BCAux::kInfiniteRange)
        return (n == 1) ? 0 : std::numeric_limits<double>::infinity();

    if (r == BCAux::kNegativeInfiniteRange)
        return ((n == 1) ? -1 : 1) * std::numeric_limits<double>::infinity();

    if (r == BCAux::kPositiveInfiniteRange)
        return std::numeric_limits<double>::infinity();

    double rm = 0;
    for (unsigned i = 0; i <= n; ++i)
        rm += pow(xmin, i) * pow(xmax, n - i);
    return rm / (n + 1);
}

// ---------------------------------------------------------
double BCConstantPrior::GetRandomValue(double xmin, double xmax, TRandom* const R)
{
    if (!R)
        return std::numeric_limits<double>::quiet_NaN();
    return xmin + R->Rndm() * xmax;
}


