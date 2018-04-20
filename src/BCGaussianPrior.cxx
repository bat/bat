/*
 * Copyright (C) 2007-2018, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

#include "BCGaussianPrior.h"
#include "BCAux.h"

#include <TMath.h>

// ---------------------------------------------------------
BCGaussianPrior::BCGaussianPrior(double mean, double sigma)
    : BCPrior()
    , fMean(mean)
    , fSigma(sigma)
{
}

// ---------------------------------------------------------
double BCGaussianPrior::GetRawMoment(unsigned n, double xmin, double xmax)
{
    if (n == 0)
        return BCPrior::GetRawMoment(n, xmin, xmax);

    BCAux::BCRange r = BCAux::RangeType(xmin, xmax);

    if (r == BCAux::kReverseRange)
        return GetRawMoment(n, xmax, xmin);

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

    double amin = (xmin - fMean) / fSigma;
    double amax = (xmax - fMean) / fSigma;
    double Z = erf(amax / sqrt(2)) - erf(amin / sqrt(2));

    if (n == 1)
        return fMean - fSigma * sqrt(2 / M_PI) * (exp(-amax * amax / 2) - exp(-amin * amin / 2)) / Z;

    // else n==2
    double aphi_min = (r == BCAux::kNegativeInfiniteRange) ? 0 : (fMean + xmin) * exp(-amin * amin / 2);
    double aphi_max = (r == BCAux::kPositiveInfiniteRange) ? 0 : (fMean + xmax) * exp(-amax * amax / 2);

    return fMean * fMean + fSigma * fSigma - fSigma * sqrt(2 / M_PI) * (aphi_max - aphi_min) / Z;
}

// ---------------------------------------------------------
double BCGaussianPrior::GetIntegral(double xmin, double xmax)
{
    switch (BCAux::RangeType(xmin, xmax)) {

        case BCAux::kReverseRange:
            return -GetIntegral(xmax, xmin);

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
            return std::numeric_limits<double>::quiet_NaN();
    }
}
