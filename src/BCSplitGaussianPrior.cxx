/*
 * Copyright (C) 2007-2015, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

#include "BCSplitGaussianPrior.h"
#include "BCAux.h"

#include <cmath>

#include <TF1.h>
#include <TMath.h>


// ---------------------------------------------------------
BCSplitGaussianPrior::BCSplitGaussianPrior(double mode, double sigma_below, double sigma_above)
    : BCPrior()
    , fMode(mode)
    , fSigmaBelow(sigma_below)
    , fSigmaAbove(sigma_above)
{
}

// ---------------------------------------------------------
BCSplitGaussianPrior::BCSplitGaussianPrior(const BCSplitGaussianPrior& other)
    : BCPrior(other)
    , fMode(other.fMode)
    , fSigmaBelow(other.fSigmaBelow)
    , fSigmaAbove(other.fSigmaAbove)
{
}

// ---------------------------------------------------------
double BCSplitGaussianPrior::GetLogPrior(double x)
{
    if (x > fMode)
        return -0.5 * (x - fMode) * (x - fMode) / fSigmaAbove / fSigmaAbove + 0.5 * log(2 / M_PI) - log(fSigmaAbove + fSigmaBelow);
    return     -0.5 * (x - fMode) * (x - fMode) / fSigmaBelow / fSigmaBelow + 0.5 * log(2 / M_PI) - log(fSigmaAbove + fSigmaBelow);
}

// ---------------------------------------------------------
double BCSplitGaussianPrior::GetRawMoment(unsigned n, double xmin, double xmax)
{
    if (n == 0 or n > 2)
        return BCPrior::GetRawMoment(n, xmin, xmax);

    BCAux::BCRange r = BCAux::RangeType(xmin, xmax);

    if (r == BCAux::kReverseRange)
        return GetRawMoment(n, xmax, xmin);

    if (r == BCAux::kEmptyRange)
        return (n == 1) ? xmin : 0;

    // if (r == BCAux::kInfiniteRange) {
    //     if (n == 1)
    //         return fMode + sqrt(2/M_PI)*(fSigmaAbove-fSigmaBelow);
    //     // else n==2
    //     return fMode * fMode + (pow(fSigmaBelow, 3) + pow(fSigmaAbove, 3) + 2 * fMode * sqrt(2 / M_PI) * (pow(fSigmaAbove, 2) - pow(fSigmaBelow, 2))) / (fSigmaBelow + fSigmaAbove);
    // }

    double smin = (xmin <= fMode) ? fSigmaBelow : fSigmaAbove;
    double smax = (xmax <= fMode) ? fSigmaBelow : fSigmaAbove;

    double erf_min = (r == BCAux::kNegativeInfiniteRange) ? -1 : TMath::Erf((xmin - fMode) / smin / sqrt(2));
    double erf_max = (r == BCAux::kPositiveInfiniteRange) ? +1 : TMath::Erf((xmax - fMode) / smax / sqrt(2));

    double phi_min = (r == BCAux::kNegativeInfiniteRange) ? 0 : exp(-0.5 * (xmin - fMode) * (xmin - fMode) / smin / smin);
    double phi_max = (r == BCAux::kPositiveInfiniteRange) ? 0 : exp(-0.5 * (xmax - fMode) * (xmax - fMode) / smax / smax);

    double I = smax * erf_max - smin * erf_min;

    // first moment
    double Ex = fMode + sqrt(2 / M_PI) * (smax * smax * (1 - phi_max) - smin * smin * (1 - phi_min)) / I;

    if (n == 1)
        return Ex;

    phi_min *= (r == BCAux::kInfiniteRange or r == BCAux::kNegativeInfiniteRange) ? 0 : sqrt(2 / M_PI) * (xmin - fMode) / smin;
    phi_max *= (r == BCAux::kInfiniteRange or r == BCAux::kPositiveInfiniteRange) ? 0 : sqrt(2 / M_PI) * (xmax - fMode) / smax;

    // second moment
    double Ex2 = 2 * fMode * Ex - fMode * fMode + (smax * smax * smax * (erf_max - phi_max) - smin * smin * smin * (erf_min - phi_min)) / I;

    return Ex2;
}

// ---------------------------------------------------------
double BCSplitGaussianPrior::GetIntegral(double xmin, double xmax)
{
    BCAux::BCRange r = BCAux::RangeType(xmin, xmax);

    if (r == BCAux::kReverseRange)
        return -GetIntegral(xmax, xmin);

    if (r == BCAux::kEmptyRange)
        return 0;

    if (r == BCAux::kInfiniteRange)
        return 1;

    double smin = (xmin <= fMode) ? fSigmaBelow : fSigmaAbove;
    double smax = (xmax <= fMode) ? fSigmaBelow : fSigmaAbove;

    double erf_min = (r == BCAux::kNegativeInfiniteRange) ? -1 : TMath::Erf((xmin - fMode) / smin / sqrt(2));
    double erf_max = (r == BCAux::kPositiveInfiniteRange) ? +1 : TMath::Erf((xmax - fMode) / smax / sqrt(2));

    return (smax * erf_max - smin * erf_min) / (fSigmaAbove + fSigmaBelow);
}

