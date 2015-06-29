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
BCTH1Prior::BCTH1Prior(TH1& h, bool interpolate)
    : BCPrior()
    , fPriorHistogram(h)
    , fInterpolate(interpolate)
{
    double integral = fPriorHistogram.Integral("width");
    if (integral != 0)
        fPriorHistogram.Scale(1. / integral);
}

// ---------------------------------------------------------
BCTH1Prior::BCTH1Prior(const BCTH1Prior& other)
    : BCPrior(other)
    , fPriorHistogram(other.fPriorHistogram)
    , fInterpolate(other.fInterpolate)
{
    double integral = fPriorHistogram.Integral("width");
    if (integral != 0)
        fPriorHistogram.Scale(1. / integral);
}

// ---------------------------------------------------------
BCTH1Prior::~BCTH1Prior()
{
}

// ---------------------------------------------------------
BCTH1Prior& BCTH1Prior::operator=(const BCTH1Prior& rhs)
{
    BCTH1Prior temp(rhs);
    swap(*this, temp);
    return *this;
}

// ---------------------------------------------------------
void swap(BCTH1Prior& A, BCTH1Prior& B)
{
    swap(static_cast<BCPrior&>(A), static_cast<BCPrior&>(B));
    TH1& temp(A.fPriorHistogram);
    A.fPriorHistogram.Copy(B.fPriorHistogram);
    B.fPriorHistogram.Copy(temp);
    std::swap(A.fInterpolate, B.fInterpolate);
}

// ---------------------------------------------------------
bool BCTH1Prior::IsValid() const
{
    if (fPriorHistogram.GetDimension() != 1)
        return false;
    double integral = 0;
    for (int i = 1; i <= fPriorHistogram.GetNbinsX(); ++i) {
        if (fPriorHistogram.GetBinContent(i) < 0)
            return false;
        integral += fPriorHistogram.GetBinContent(i);
    }
    if (integral == 0)
        return false;
    return true;
}

// ---------------------------------------------------------
double BCTH1Prior::GetLogPrior(double x) const
{
    double p = GetPrior(x, false);
    if (p > 0)
        return log(p);
    if (p == 0)
        return -std::numeric_limits<double>::infinity();
    return std::numeric_limits<double>::quiet_NaN();
}

// ---------------------------------------------------------
double BCTH1Prior::GetMode(double /*xmin*/, double /*xmax*/) const
{
    return fPriorHistogram.GetXaxis()->GetBinCenter(fPriorHistogram.GetMaximumBin());
}

// ---------------------------------------------------------
double BCTH1Prior::GetRawMoment(unsigned n, double xmin, double xmax) const
{
    if (n == 0)
        return 0;
    if (n == 1)
        return fPriorHistogram.GetMean();
    if (n == 2)
        return fPriorHistogram.GetMean() * fPriorHistogram.GetMean() + fPriorHistogram.GetRMS() * fPriorHistogram.GetRMS();
    return BCPrior::GetRawMoment(n, xmin, xmax);
}

// ---------------------------------------------------------
double BCTH1Prior::GetStandardizedMoment(unsigned n, double xmin, double xmax) const
{
    if (n == 0)
        return 0;
    if (n == 1)
        return 1;
    if (n == 2)
        return GetSkewness(xmin, xmax);
    if (n == 3)
        return GetKurtosis(xmin, xmax);
    return BCPrior::GetStandardizedMoment(n, xmin, xmax);
}

// ---------------------------------------------------------
double BCTH1Prior::GetIntegral(double xmin, double xmax) const
{
    xmin = std::max(xmin, fPriorHistogram.GetXaxis()->GetXmin());
    xmax = std::min(xmax, fPriorHistogram.GetXaxis()->GetXmax());
    int bmin = fPriorHistogram.FindFixBin(xmin);
    int bmax = fPriorHistogram.FindFixBin(xmax);
    double I = fPriorHistogram.Integral(bmin, bmax, "width");
    I -= fPriorHistogram.GetBinContent(xmin) * (xmin - fPriorHistogram.GetXaxis()->GetBinLowEdge(bmin));
    I -= fPriorHistogram.GetBinContent(xmax) * (fPriorHistogram.GetXaxis()->GetBinUpEdge(bmax) - xmax);
    return I;
}
