/*
 * Copyright (C) 2007-2018, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

// ---------------------------------------------------------

#include "BCAux.h"
#include "BCLog.h"
#include "BCPrior.h"
#include "BCConstantPrior.h"

#include <TMath.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TRandom.h>

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
double BCPrior::GetPrior(double x, bool normalize)
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
double BCPrior::GetMode(double xmin, double xmax)
{
    BCAux::MakeFinite(xmin, xmax);
    return fPriorFunction.GetMaximumX(xmin, xmax);
}

// ---------------------------------------------------------
double BCPrior::GetRawMoment(unsigned n, double xmin, double xmax)
{
    if (n == 0)
        return 1;
    BCAux::MakeFinite(xmin, xmax);
    return fPriorFunction.Moment(static_cast<double>(n), xmin, xmax);
}

// ---------------------------------------------------------
double BCPrior::GetIntegral(double xmin, double xmax)
{
    BCAux::MakeFinite(xmin, xmax);
    return fPriorFunction.Integral(xmin, xmax);
}

// ---------------------------------------------------------
double BCPrior::GetCentralMoment(unsigned n, double xmin, double xmax)
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
double BCPrior::GetStandardizedMoment(unsigned n, double xmin, double xmax)
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
double BCPrior::GetRandomValue(double xmin, double xmax, TRandom* const R)
{
    (void)R;
    return fPriorFunction.GetRandom(xmin, xmax);
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
BCH1D BCPrior::GetBCH1D(TH1* bins, const std::string& name)
{
    BCH1D bch1;
    if (!IsValid())
        return bch1;

    TH1* h = BCAux::OwnClone(bins, name);

    h->Add(&GetFunction(), 1, "I");

    bch1 = h;
    bch1.SetLocalMode(GetMode(h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax()));

    return bch1;
}

// ---------------------------------------------------------
BCH2D BCPrior::GetBCH2D(BCPrior* ordinate, TH2* bins, const std::string& name)
{
    BCH2D bch2;

    if (!ordinate || !ordinate->IsValid())
        return bch2;

    BCH1D bch_x = GetBCH1D(bins->ProjectionX(), "tempx");
    if (!bch_x.Valid())
        return bch2;

    BCH1D bch_y = ordinate->GetBCH1D(bins->ProjectionY(), "tempy");
    if  (!bch_y.Valid())
        return bch2;

    TH2* h = BCAux::OwnClone(bins, name);

    // get x bins:
    std::vector<double> bin_edges_x(bch_x.GetHistogram()->GetNbinsX() + 1, 0);
    bch_x.GetHistogram()->GetXaxis()->GetLowEdge(&bin_edges_x[0]);
    bin_edges_x.back() = bch_x.GetHistogram()->GetXaxis()->GetXmax();
    // get y bins
    std::vector<double> bin_edges_y(bch_y.GetHistogram()->GetNbinsX() + 1, 0);
    bch_y.GetHistogram()->GetXaxis()->GetLowEdge(&bin_edges_y[0]);
    bin_edges_y.back() = bch_y.GetHistogram()->GetXaxis()->GetXmax();
    // set into bins histogram
    bins->SetBins(bin_edges_x.size() - 1, &bin_edges_x[0], bin_edges_y.size() - 1, &bin_edges_y[0]);

    // set content
    for (int x = 1; x <= bch_x.GetHistogram()->GetNbinsX(); ++x)
        for (int y = 1; y <= bch_y.GetHistogram()->GetNbinsX(); ++y)
            h->SetBinContent(x, y, bch_x.GetHistogram()->GetBinContent(x) * bch_x.GetHistogram()->GetBinContent(y));

    bch2 = h;

    return bch2;
}
