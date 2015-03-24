/*
 * Copyright (C) 2007-2015, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

// ---------------------------------------------------------

#include "BCPrior.h"
#include "BCLog.h"
#include "BCAux.h"

#include <TMath.h>
#include <TF1.h>
#include <TRandom.h>

#include <iostream>

// ---------------------------------------------------------
BCPrior::BCPrior()
    : fPriorFunction(new TF1("prior_interal_f1", this, &BCPrior::GetPriorForROOT, 0, 0, 1)) //-std::numeric_limits<double>::max(),std::numeric_limits<double>::max()))
    , fIntegral(1)
{
}

// ---------------------------------------------------------
BCPrior::BCPrior(const TF1* const f)
    : fPriorFunction(NULL)
    , fIntegral(1)
{
    if (f)
        fPriorFunction = new TF1(*f);
}

// ---------------------------------------------------------
BCPrior::BCPrior(const BCPrior& other)
    : fPriorFunction(new TF1("prior_interal_f1", this, &BCPrior::GetPriorForROOT, 0, 0, 1)) //-std::numeric_limits<double>::max(),std::numeric_limits<double>::max()))
    , fIntegral(other.fIntegral)
{
}

// ---------------------------------------------------------
BCPrior::~BCPrior()
{
    delete fPriorFunction;
}

// ---------------------------------------------------------
double BCPrior::GetPrior(double x) const
{
    double lp = GetLogPrior(x);
    if (std::isfinite(lp))
        return exp(lp);
    if (lp < 0)
        return 0;										// exp(-inf)
    return std::numeric_limits<double>::infinity(); // exp(inf)
}

// ---------------------------------------------------------
double BCPrior::GetMode(double xmin, double xmax) const
{
    if (!fPriorFunction)
        return std::numeric_limits<double>::quiet_NaN();
    BCAux::MakeFinite(xmin, xmax);
    return fPriorFunction->GetMaximumX(xmin, xmax);
}

// ---------------------------------------------------------
double BCPrior::GetRawMoment(unsigned n, double xmin, double xmax) const
{
    if (n == 0)
        return 1;
    if (!fPriorFunction)
        return std::numeric_limits<double>::infinity();
    BCAux::MakeFinite(xmin, xmax);
    return fPriorFunction->Moment(static_cast<double>(n), xmin, xmax);
}

// ---------------------------------------------------------
double BCPrior::GetIntegral(double xmin, double xmax) const
{
    if (!fPriorFunction)
        return 0;
    BCAux::MakeFinite(xmin, xmax);
    return fPriorFunction->Integral(xmin, xmax);
}

// ---------------------------------------------------------
double BCPrior::GetCentralMoment(unsigned n, double xmin, double xmax) const
{
    if (n == 0)
        return std::numeric_limits<double>::infinity();

    if (n == 1)
        return 0;

    double mean = GetRawMoment(1, xmin, xmax);
    if (!std::isfinite(mean))
        return std::numeric_limits<double>::infinity();

    double cm = 0;
    for (unsigned i = n; i > 1; --i) {
        double rm = GetRawMoment(i, xmin, xmax);
        if (!std::isfinite(rm))
            return std::numeric_limits<double>::infinity();
        cm += TMath::Binomial(n, i) * rm * pow(-mean, n - i);
    }
    cm -= (n - 1) * pow(-mean, n);
    return cm;
}

// ---------------------------------------------------------
double BCPrior::GetStandardizedMoment(unsigned n, double xmin, double xmax) const
{
    double variance = GetVariance(xmin, xmax);
    if (!std::isfinite(variance))
        return std::numeric_limits<double>::infinity();

    double cm = GetCentralMoment(n, xmin, xmax);
    if (!std::isfinite(cm))
        return std::numeric_limits<double>::infinity();

    return cm / pow(variance, n / 2.);
}

// ---------------------------------------------------------
double BCPrior::GetRandomValue(double xmin, double xmax, TRandom* const R) const
{
    (void)R;
    if (!fPriorFunction)
        return std::numeric_limits<double>::quiet_NaN();
    return fPriorFunction->GetRandom(xmin, xmax);
}

// ---------------------------------------------------------
double BCPrior::GetRandomValueGaussian(double xmin, double xmax, TRandom* const R, double expansion_factor, unsigned N, bool over_range) const
{
    if (!R) {
        BCLog::OutError("BCPrior::GetRandomValueGaussian : no random number generator provided.");
        return std::numeric_limits<double>::quiet_NaN();
    }

    double m = std::numeric_limits<double>::quiet_NaN();
    double s = std::numeric_limits<double>::quiet_NaN();

    BCAux::BCRange r = BCAux::RangeType(xmin, xmax);

    if (!over_range) {
        m = GetMean();
        s = GetStandardDeviation();
        // if could not calculate mean in infinite range
        if (!std::isfinite(m) and r != BCAux::kInfiniteRange) {
            BCLog::OutWarning("BCPrior::GetRandomValueGaussian : Infinite range failed, using confined range!");
            over_range = true;
        }
    }

    if (over_range) {
        m = GetMean(xmin, xmax);
        s = GetStandardDeviation(xmin, xmax);
    }

    // if could not calculate mean and sigma
    if (!std::isfinite(m)) {
        BCLog::OutError("BCPrior::GetRandomValueGaussian : could not calculate prior's mean!");
        return std::numeric_limits<double>::quiet_NaN();
    }

    // increase standard deviation by expansion_factor
    s *= expansion_factor;

    // if standard deviation is infinite, use uniform distribution
    if (!std::isfinite(s))
        return (r == BCAux::kFiniteRange) ? xmin + R->Rndm() * (xmax - xmin) : std::numeric_limits<double>::quiet_NaN();

    // if standard deviation is zero, use delta(x-mean)
    if (s == 0)
        return (m >= xmin and m <= xmax) ? m : std::numeric_limits<double>::quiet_NaN();

    // check overlap of range with gaussian and issue warning if necessary
    if (m<xmin or m>xmax) {
        double d = ((m < xmin) ? xmin - m : m - xmax) / s;
        if (d > 4)
            BCLog::OutWarning("BCrior::GetRandomValueGaussian : requested range beyond 3 sigma of prior mean, this may take a while...");
        else if ((xmax - xmin) / s < 1e-4)
            BCLog::OutWarning("BCPrior::GetRandomValueGaussian : requested range very small compared to prior width, this may take a while...");
    }

    // repeat until value inside range is selected
    // or maximum number of tries is exhausted
    for (unsigned n = 0; n < N; ++n) {
        double x = R->Gaus(m, s);
        if (x >= xmin and x <= xmax)
            return x;
        x = 2 * m - s;
        if (x >= xmin and x <= xmax)
            return x;
    }
    return std::numeric_limits<double>::quiet_NaN();
}
