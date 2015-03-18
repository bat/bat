/*
 * Copyright (C) 2007-2015, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

#include "BCGaussianPrior.h"
#include "BCAux.h"

#include <TMath.h>

#include <iostream>

// ---------------------------------------------------------
BCGaussianPrior::BCGaussianPrior(double mean, double sigma)
    : BCPrior()
    , fMean(mean)
    , fSigma(sigma)
{
}

// ---------------------------------------------------------
BCGaussianPrior::BCGaussianPrior(const BCGaussianPrior& other)
    : BCPrior(other)
    , fMean(other.fMean)
    , fSigma(other.fSigma)
{}

// // ---------------------------------------------------------
BCGaussianPrior::~BCGaussianPrior()
{
}

// ---------------------------------------------------------
double BCGaussianPrior::GetRawMoment(unsigned n, double xmin, double xmax) const
{
    if (n == 0)
        return BCPrior::GetRawMoment(n, xmin, xmax);

    BCAux::BCRange r = BCAux::RangeType(xmin, xmax);

    if (r == BCAux::kInfiniteRange) {
        if (n == 1)
            return fMean;
        if (n == 2)
            return fMean * fMean + fSigma * fSigma;
        if (n == 3)
            return fMean * (fMean * fMean + 3 * fSigma * fSigma);
        if (n == 4)
            return pow(fMean, 4) + 6 * fMean * fMean * fSigma * fSigma + 3 * pow(fSigma, 4);
        return BCPrior::GetRawMoment(n, xmin, xmax);
    }

    if (r == BCAux::kEmptyRange)
        return (n == 1) ? xmin : 0;

    if (n > 2)
        return BCPrior::GetRawMoment(n, xmin, xmax);

    double gmin = (r == BCAux::kNegativeInfiniteRange) ? 0 : exp(-0.5 * (xmin - fMean) * (xmin - fMean) / fSigma / fSigma) / fSigma / sqrt(2 * M_PI);
    double gmax = (r == BCAux::kPositiveInfiniteRange) ? 0 : exp(-0.5 * (xmax - fMean) * (xmax - fMean) / fSigma / fSigma) / fSigma / sqrt(2 * M_PI);

    if (n == 1)
        return fMean - fSigma * fSigma * (gmax - gmin) / GetIntegral(xmin, xmax) / 2;

    // else n==2
    gmin *= (r == BCAux::kNegativeInfiniteRange) ? 0 : xmin;
    gmax *= (r == BCAux::kPositiveInfiniteRange) ? 0 : xmax;
    return fMean * GetMean(xmin, xmax) + fSigma * fSigma * (1 - (gmax - gmin) / GetIntegral(xmin, xmax));
}

// ---------------------------------------------------------
double BCGaussianPrior::GetIntegral(double xmin, double xmax) const
{
    switch (BCAux::RangeType(xmin, xmax)) {

        case BCAux::kFiniteRange:
            return (TMath::Erf((xmax - fMean) / fSigma / sqrt(2)) - TMath::Erf((xmin - fMean) / fSigma / sqrt(2))) / 2;

        case BCAux::kNegativeInfiniteRange:
            return (1 + TMath::Erf((xmax - fMean) / fSigma / sqrt(2))) / 2;

        case BCAux::kPositiveInfiniteRange:
            return (1 - TMath::Erf((xmin - fMean) / fSigma / sqrt(2))) / 2;

        case BCAux::kInfiniteRange:
            return 1;

        case BCAux::kEmptyRange:
            return 0;

        default:
            return std::numeric_limits<double>::infinity();
    }
}
