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
BCSplitGaussianPrior::BCSplitGaussianPrior(double mean, double sigma_below, double sigma_above)
    : BCPrior()
    , fMean(mean)
    , fSigmaBelow(sigma_below)
    , fSigmaAbove(sigma_above)
{
}

// ---------------------------------------------------------
BCSplitGaussianPrior::BCSplitGaussianPrior(const BCSplitGaussianPrior& other)
    : BCPrior(other)
    , fMean(other.fMean)
    , fSigmaBelow(other.fSigmaBelow)
    , fSigmaAbove(other.fSigmaAbove)
{
}

// ---------------------------------------------------------
BCSplitGaussianPrior::~BCSplitGaussianPrior()
{
}

// ---------------------------------------------------------
double BCSplitGaussianPrior::GetLogPrior(double x) const
{
    if (x > fMean)
        return -0.5*(x-fMean)*(x-fMean)/fSigmaAbove/fSigmaAbove + 0.5*log(2/M_PI) - log(fSigmaAbove + fSigmaBelow);
    return     -0.5*(x-fMean)*(x-fMean)/fSigmaBelow/fSigmaBelow + 0.5*log(2/M_PI) - log(fSigmaAbove + fSigmaBelow);
}

// ---------------------------------------------------------
double BCSplitGaussianPrior::GetRawMoment(unsigned n, double xmin, double xmax) const
{
    if (n == 0 or n > 2)
        return BCPrior::GetRawMoment(n, xmin, xmax);

    BCAux::BCRange r = BCAux::RangeType(xmin, xmax);

    if (r == BCAux::kReverseRange)
        return GetRawMoment(n,xmax,xmin);

    if (r == BCAux::kEmptyRange)
        return (n == 1) ? xmin : 0;

    if (r == BCAux::kInfiniteRange) {
        if (n == 1)
            return fMean + sqrt(2/M_PI)*(fSigmaAbove-fSigmaBelow);
        // else n==2
        return fMean * fMean + (pow(fSigmaBelow, 3) + pow(fSigmaAbove, 3) + 2 * fMean * sqrt(2 / M_PI) * (pow(fSigmaAbove, 2) - pow(fSigmaBelow, 2))) / (fSigmaBelow + fSigmaAbove);
    }

    double smin = (xmin<=fMean) ? fSigmaBelow : fSigmaAbove;
    double smax = (xmax<=fMean) ? fSigmaBelow : fSigmaAbove;

    double erf_min = (r == BCAux::kNegativeInfiniteRange) ? -1 : TMath::Erf((xmin-fMean)/smin/sqrt(2));
    double erf_max = (r == BCAux::kPositiveInfiniteRange) ? +1 : TMath::Erf((xmax-fMean)/smax/sqrt(2));

    double phi_min = (r == BCAux::kNegativeInfiniteRange) ? 0 : exp(-0.5*(xmin-fMean)*(xmin-fMean)/smin/smin);
    double phi_max = (r == BCAux::kPositiveInfiniteRange) ? 0 : exp(-0.5*(xmax-fMean)*(xmax-fMean)/smax/smax);

    if (n == 1)
        return fMean + sqrt(2/M_PI)*(smax*smax*(1-phi_max)-smin*smin*(1-phi_min))/(smax*erf_max-smin*erf_min);

    phi_min *= (r == BCAux::kNegativeInfiniteRange) ? 0 : (fMean+xmin);
    phi_max *= (r == BCAux::kPositiveInfiniteRange) ? 0 : (fMean+xmax);

    // else n == 2
    return fMean*fMean + sqrt(2/M_PI)*(smax*smax*(erf_max+2*fMean-phi_max)-smin*smin*(erf_min+2*fMean-phi_min))/(smax*erf_max-smin*erf_min);
    
    // double phi_min = (r == BCAux::kNegativeInfiniteRange) ? 0 : exp(-0.5*(xmin - fMean) / fSigmaBelow, 2));
    // double phi_max = (r == BCAux::kPositiveInfiniteRange) ? 0 : exp(-0.5*pow((xmax - fMean) / fSigmaAbove, 2));

    // double erf_min = (r == BCAux::kNegativeInfiniteRange) ? -1 : TMath::Erf((xmin - fMean) / fSigmaBelow / sqrt(2));
    // double erf_max = (r == BCAux::kPositiveInfiniteRange) ? +1 : TMath::Erf((xmax - fMean) / fSigmaAbove / sqrt(2));

    // if (n == 1)
    //     return fMean + sqrt(2 / M_PI) * (pow(fSigmaAbove, 2) * (1 - phi_max) - pow(fSigmaBelow, 2) * (1 - phi_min)) / (fSigmaAbove * erf_max - fSigmaBelow * erf_min);

    // else n==2
    // phi_min *= (r == BCAux::kNegativeInfiniteRange) ? 1 : (1 + xmin / fMean);
    // phi_min *= (r == BCAux::kPositiveInfiniteRange) ? 1 : (1 + xmax / fMean);
    // return fMean * fMean + (sqrt(2 / M_PI) * fMean * (pow(fSigmaAbove, 2) * (2 - phi_max) - pow(fSigmaBelow, 2) * (2 - phi_min)) + pow(fSigmaAbove, 3) * erf_max - pow(fSigmaBelow, 3) * erf_min) / (fSigmaAbove * erf_max - fSigmaBelow * erf_min);
}

// ---------------------------------------------------------
double BCSplitGaussianPrior::GetIntegral(double xmin, double xmax) const
{
    BCAux::BCRange r = BCAux::RangeType(xmin, xmax);

    if (r == BCAux::kReverseRange)
        return -GetIntegral(xmax,xmin);

    if (r == BCAux::kEmptyRange)
        return 0;

    if (r == BCAux::kInfiniteRange)
        return 1;
    
    double smin = (xmin<=fMean) ? fSigmaBelow : fSigmaAbove;
    double smax = (xmax<=fMean) ? fSigmaBelow : fSigmaAbove;

    double erf_min = (r == BCAux::kNegativeInfiniteRange) ? -1 : TMath::Erf((xmin-fMean)/smin/sqrt(2));
    double erf_max = (r == BCAux::kPositiveInfiniteRange) ? +1 : TMath::Erf((xmax-fMean)/smax/sqrt(2));
    
    return (smax*erf_max - smin*erf_min)/(fSigmaAbove+fSigmaBelow);
}

