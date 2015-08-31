/*
 * Copyright (C) 2007-2015, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

#include <config.h>

#include "BCPriorTH1.h"

#include <cmath>

// ---------------------------------------------------------
BCPriorTH1::BCPriorTH1(TH1& h, bool interpolate)
    : BCPrior(),
      fPriorHistogram(h),
      fInterpolate(interpolate)
{
    NormalizeHistogram();
}

// ---------------------------------------------------------
BCPriorTH1::BCPriorTH1(TH1* h, bool interpolate)
    : BCPrior(),
      fPriorHistogram(*h),
      fInterpolate(interpolate)
{
    NormalizeHistogram();
}

// ---------------------------------------------------------
BCPriorTH1::BCPriorTH1(const BCPriorTH1& other)
    : BCPrior(other),
      fPriorHistogram(other.fPriorHistogram),
      fInterpolate(other.fInterpolate)
{
    NormalizeHistogram();
}

// ---------------------------------------------------------
BCPriorTH1& BCPriorTH1::operator=(BCPriorTH1 rhs)
{
    swap(*this, rhs);
    return *this;
}

// ---------------------------------------------------------
void swap(BCPriorTH1& A, BCPriorTH1& B)
{
    swap(static_cast<BCPrior&>(A), static_cast<BCPrior&>(B));
    TH1& temp(A.fPriorHistogram);
#if ROOTVERSION < 5034025
    // ROOT versions before 5.34/25 don't have a public copy function
    new (&A) TH1(B.fPriorHistogram);
    new (&B) TH1(temp);
#else
    // newer versions do
    A.fPriorHistogram.Copy(B.fPriorHistogram);
    B.fPriorHistogram.Copy(temp);
#endif
    std::swap(A.fInterpolate, B.fInterpolate);
}

// ---------------------------------------------------------
bool BCPriorTH1::IsValid() const
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
void BCPriorTH1::NormalizeHistogram()
{
    double integral = 0;
    if (fInterpolate)
        integral = GetFunction().Integral(fPriorHistogram.GetXaxis()->GetXmin(), fPriorHistogram.GetXaxis()->GetXmax());
    else
        integral = fPriorHistogram.Integral("width");

    if (integral != 0)
        fPriorHistogram.Scale(1. / integral);
}

// ---------------------------------------------------------
double BCPriorTH1::GetLogPrior(double x)
{
    double p = GetPrior(x, false);
    if (p > 0)
        return log(p);
    if (p == 0)
        return -std::numeric_limits<double>::infinity();
    return std::numeric_limits<double>::quiet_NaN();
}

// ---------------------------------------------------------
double BCPriorTH1::GetMode(double /*xmin*/, double /*xmax*/)
{
    return fPriorHistogram.GetXaxis()->GetBinCenter(fPriorHistogram.GetMaximumBin());
}

// ---------------------------------------------------------
double BCPriorTH1::GetRawMoment(unsigned n, double xmin, double xmax)
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
double BCPriorTH1::GetStandardizedMoment(unsigned n, double xmin, double xmax)
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
double BCPriorTH1::GetIntegral(double xmin, double xmax)
{
    xmin = std::max(xmin, fPriorHistogram.GetXaxis()->GetXmin());
    xmax = std::min(xmax, fPriorHistogram.GetXaxis()->GetXmax());

    // if interpolating, use built in function
    if (fInterpolate)
        return const_cast<TF1*>(&GetFunction())->Integral(xmin, xmax);

    // else calculate directly from histogram
    int bmin = fPriorHistogram.FindFixBin(xmin);
    int bmax = fPriorHistogram.FindFixBin(xmax);
    double I = fPriorHistogram.Integral(bmin, bmax, "width");
    I -= fPriorHistogram.GetBinContent(xmin) * (xmin - fPriorHistogram.GetXaxis()->GetBinLowEdge(bmin));
    I -= fPriorHistogram.GetBinContent(xmax) * (fPriorHistogram.GetXaxis()->GetBinUpEdge(bmax) - xmax);
    return I;
}

// ---------------------------------------------------------
BCH1D BCPriorTH1::GetBCH1D(TH1* bins, const char* name)
{
    // if not interpolating, use actual histogram binning
    if (!fInterpolate) {
        std::vector<double> bin_edges(fPriorHistogram.GetNbinsX() + 1, 0);
        fPriorHistogram.GetXaxis()->GetLowEdge(&bin_edges[0]);
        bin_edges.back() = fPriorHistogram.GetXaxis()->GetXmax();
        bins->SetBins(bin_edges.size() - 1, &bin_edges[0]);
    }
    return BCPrior::GetBCH1D(bins, name);
}
