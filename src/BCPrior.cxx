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
#include <TH1.h>
#include <TRandom.h>

#include <iostream>

// ---------------------------------------------------------
BCPrior::BCPrior()
    : fPriorFunction("prior_interal_f1", this, &BCPrior::GetPriorForROOT, 0, 0, 1) //-std::numeric_limits<double>::max(),std::numeric_limits<double>::max()))
    , fLogIntegral(0)
{
}

// ---------------------------------------------------------
BCPrior::BCPrior(const BCPrior& other)
    : fPriorFunction("prior_interal_f1", this, &BCPrior::GetPriorForROOT, 0, 0, 1) //-std::numeric_limits<double>::max(),std::numeric_limits<double>::max()))
    , fLogIntegral(other.fLogIntegral)
{
}

// ---------------------------------------------------------
BCPrior::~BCPrior()
{
}

// ---------------------------------------------------------
void swap(BCPrior& A, BCPrior& B)
{
    TF1 temp = A.fPriorFunction;
    A.fPriorFunction = B.fPriorFunction;
    B.fPriorFunction = temp;
    std::swap(A.fLogIntegral, B.fLogIntegral);
}

// ---------------------------------------------------------
double BCPrior::GetPrior(double x, bool normalize) const
{
    double lp = GetLogPrior(x);
    if (normalize)
        lp -= fLogIntegral;

    if (std::isfinite(lp))
        return exp(lp);
    if (lp < 0)
        return 0;										// exp(-inf)
    return std::numeric_limits<double>::infinity(); // exp(inf)
}

// ---------------------------------------------------------
double BCPrior::GetMode(double xmin, double xmax) const
{
    BCAux::MakeFinite(xmin, xmax);
    return fPriorFunction.GetMaximumX(xmin, xmax);
}

// ---------------------------------------------------------
double BCPrior::GetRawMoment(unsigned n, double xmin, double xmax) const
{
    if (n == 0)
        return 1;
    BCAux::MakeFinite(xmin, xmax);
    return const_cast<TF1*>(&fPriorFunction)->Moment(static_cast<double>(n), xmin, xmax);
}

// ---------------------------------------------------------
double BCPrior::GetIntegral(double xmin, double xmax) const
{
    BCAux::MakeFinite(xmin, xmax);
    return const_cast<TF1*>(&fPriorFunction)->Integral(xmin, xmax);
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
    return const_cast<TF1*>(&fPriorFunction)->GetRandom(xmin, xmax);
}

// ---------------------------------------------------------
void BCPrior::SetFunctionRange(double xmin, double xmax)
{
    fPriorFunction.SetRange(xmin, xmax);
    CalculateAndStoreIntegral(xmin, xmax);
}

// ---------------------------------------------------------
void BCPrior::FillHistogramByCenterValue(TH1* h)
{
    if (!h)
        return;
    for (int i = 1; i <= h->GetNbinsX(); ++i)
        h -> SetBinContent(i, GetPrior(h->GetXaxis()->GetBinCenter(i)));
}

// ---------------------------------------------------------
void BCPrior::FillHistogramByIntegral(TH1* h)
{
    if (!h)
        return;
    for (int i = 1; i <= h->GetNbinsX(); ++i)
        if (h->GetXaxis()->GetBinWidth(i) > 0)
            h -> SetBinContent(i, GetIntegral(h->GetXaxis()->GetBinLowEdge(i), h->GetXaxis()->GetBinUpEdge(i)) / h->GetXaxis()->GetBinWidth(i));
}

// ---------------------------------------------------------
BCH1D BCPrior::GetBCH1D(TH1* bins, const char* name)
{
    TH1* h = (TH1*) bins->Clone(name);
    h->Add(&GetFunction(), 1, "I");

    BCH1D bch = h;
    bch.SetLocalMode(GetMode(h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax()));

    return bch;
}
