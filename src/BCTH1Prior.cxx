/*
 * Copyright (C) 2007-2015, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

#include "BCTH1Prior.h"

#include <cmath>

// ---------------------------------------------------------
BCTH1Prior::BCTH1Prior(TH1 const* const h, bool interpolate)
    : BCPrior()
    , fPriorHistogram(NULL)
    , fInterpolate(interpolate)
{
    if (h and h->GetDimension() == 1)
        fPriorHistogram = (TH1*) h->Clone();
    double integral = fPriorHistogram->Integral("width");
    if (integral != 0)
        fPriorHistogram->Scale(1. / integral);
}

// ---------------------------------------------------------
BCTH1Prior::BCTH1Prior(const BCTH1Prior& other)
    : BCPrior(other)
    , fPriorHistogram(NULL)
    , fInterpolate(other.fInterpolate)
{
    if (other.fPriorHistogram)
        fPriorHistogram = (TH1*) other.fPriorHistogram->Clone();
    double integral = fPriorHistogram->Integral("width");
    if (integral != 0)
        fPriorHistogram->Scale(1. / integral);
}

// ---------------------------------------------------------
BCTH1Prior::~BCTH1Prior()
{
    delete fPriorHistogram;
}

// ---------------------------------------------------------
double BCTH1Prior::GetMode(double /*xmin*/, double /*xmax*/) const
{
    if (!fPriorHistogram)
        return std::numeric_limits<double>::quiet_NaN();
    return fPriorHistogram->GetXaxis()->GetBinCenter(fPriorHistogram->GetMaximumBin());
}


// ---------------------------------------------------------
double BCTH1Prior::GetRawMoment(unsigned n, double xmin, double xmax) const
{
    if (n == 0)
        return 0;
    if (n == 1 and fPriorHistogram)
        return fPriorHistogram->GetMean();
    if (n == 2 and fPriorHistogram)
        return fPriorHistogram->GetMean() * fPriorHistogram->GetMean() + fPriorHistogram->GetRMS() * fPriorHistogram->GetRMS();
    return BCPrior::GetRawMoment(n, xmin, xmax);
}

// ---------------------------------------------------------
double BCTH1Prior::GetStandardisedMoment(unsigned n, double xmin, double xmax) const
{
    if (n == 0)
        return 0;
    if (n == 1)
        return 1;
    if (n == 2 and fPriorHistogram)
        return GetSkewness(xmin, xmax);
    if (n == 3 and fPriorHistogram)
        return GetKurtosis(xmin, xmax);
    return BCPrior::GetStandardisedMoment(n, xmin, xmax);
}

// ---------------------------------------------------------
double BCTH1Prior::GetRandomValue(double /*xmin*/, double /*xmax*/, TRandom* const /*R*/) const
{
    if (!fPriorHistogram)
        return std::numeric_limits<double>::quiet_NaN();
    return fPriorHistogram->GetRandom();
}
