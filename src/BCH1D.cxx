/*
 * Copyright (C) 2007-2018, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

// ---------------------------------------------------------

#include "BCH1D.h"

#include "BCAux.h"
#include "BCLog.h"

#include <TArrow.h>
#include <TAxis.h>
#include <TH1D.h>
#include <TLegendEntry.h>
#include <TLine.h>
#include <TMarker.h>
#include <TPad.h>
#include <TROOT.h>
#include <TString.h>
#include <TStyle.h>

#include <algorithm>
#include <limits>
#include <math.h>

// ---------------------------------------------------------
BCH1D::BCH1D(const TH1* const hist)
    : BCHistogramBase(hist, 1)
    , fBandType(kSmallestInterval)
    , fNQuantiles(0)
    , fQuantileLineColor(kBlack)
    , fDrawMedian(false)
    , fDrawCentral68(true)
{
}

// ---------------------------------------------------------
BCH1D::BCH1D(const BCH1D& other)
    : BCHistogramBase(other)
{
    CopyOptions(other);
}

// ---------------------------------------------------------
void BCH1D::CopyOptions(const BCH1D& other)
{
    BCHistogramBase::CopyOptions(other);
    fBandType = other.fBandType;
    fNQuantiles = other.fNQuantiles;
    fQuantileLineColor = other.fQuantileLineColor;
    fDrawMedian = other.fDrawMedian;
    fDrawCentral68 = other.fDrawCentral68;
}

// ---------------------------------------------------------
void BCH1D::SetColorScheme(BCHColorScheme scheme)
{
    BCHistogramBase::SetColorScheme(scheme);
    SetQuantileLineColor(GetMarkerColor());
}

// ---------------------------------------------------------
double BCH1D::GetQuantile(double probability)
{
    int nquantiles = 1;
    double quantiles[1];
    double probsum[1];

    probsum[0] = probability;

    // use ROOT function to calculate quantile.
    fHistogram->GetQuantiles(nquantiles, quantiles, probsum);

    return quantiles[0];
}

// ---------------------------------------------------------
void BCH1D::CheckIntervals(std::vector<double>& intervals)
{
    // check & sort interval values
    if (fBandType != kNoBands)
        BCHistogramBase::CheckIntervals(fIntervals, (fBandType == kLowerLimit) ? -1 : +1);

    // check number of intervals values if user-specified
    if (fBandType == kUserSpecified && intervals.size() == 1) {
        BCLog::OutWarning("BCH1D::CheckIntervals : at least two intervals values must be specified for user-specified intervals. No bands will be drawn.");
        intervals.clear();
    }
}

// ---------------------------------------------------------
std::vector<double> BCH1D::DefaultIntervals(int nbands)
{
    if (nbands < 0)
        nbands = fNBands;

    std::vector<double> intervals;

    switch (fBandType) {

        case kNoBands:
        case kUserSpecified:
            return intervals;

        case kUpperLimit:
            if (nbands > 0) intervals.push_back(0.90);
            if (nbands > 1) intervals.push_back(0.95);
            if (nbands > 2) intervals.push_back(0.99);
            for (int i = 3; i < nbands; ++i) // not recommended to go (far) beyond 3
                intervals.push_back(0.9 + intervals.back() / 10.);
            return intervals;

        case kLowerLimit:
            if (nbands > 0) intervals.push_back(0.10);
            if (nbands > 1) intervals.push_back(0.05);
            if (nbands > 2) intervals.push_back(0.01);
            for (int i = 3; i < nbands; ++i) // not recommended to go (far) beyond 3
                intervals.push_back(intervals.back() / 10.);
            return intervals;

        case kCentralInterval:
        case kSmallestInterval:
        default:
            return BCHistogramBase::DefaultIntervals(nbands);

    }
}

// ---------------------------------------------------------
void BCH1D::DrawBands(const std::string& options)
{
    GetHistogram()->Draw(options.data());

    if ( fBandType == kNoBands || GetHistogram()->Integral() <= 0 )
        return;

    std::vector<double> intervals = fIntervals;
    CheckIntervals(intervals);

    if (intervals.empty())
        return;

    unsigned nbands = intervals.size() - ((fBandType == kUserSpecified) ? 1 : 0);

    // calculate bounds for smallest intervals, if using
    std::vector<std::pair<double, double> > bounds;
    if (fBandType == kSmallestInterval) {
        bounds = GetSmallestIntervalBounds(intervals, fBandOvercoverage);
        nbands = bounds.size();
    }

    // make sure enough colors have been designated
    while (nbands > fBandColors.size())
        fBandColors.push_back(fBandColors.back() - 1);

    int i0 = fROOTObjects.size();
    for (int i = nbands - 1; i >= 0; --i) {

        TH1* hist_band = NULL;
        std::string legend_text;

        if (fBandType == kSmallestInterval) {
            hist_band = BCAux::OwnClone(GetHistogram(), std::string(GetHistogram()->GetName()) + "_band");

            for (int b = 1; b <= hist_band->GetNbinsX(); ++b)
                if (hist_band->GetBinContent(b) < bounds[i].first)
                    hist_band->SetBinContent(b, 0);
            legend_text = "smallest %.1f%% interval(s)";
        } else {
            double p[2] = {0, 1};
            switch (fBandType) {

                case kUpperLimit:
                    p[1] = intervals[i];
                    legend_text = "%.0f%% upper limit";
                    break;

                case kLowerLimit:
                    p[0] = intervals[i];
                    legend_text = "%.0f%% lower limit";
                    break;

                case kUserSpecified:
                    p[0] = intervals[i];
                    p[1] = intervals[i + 1];
                    legend_text = Form("%%.1f%%%% interval from %.1f%%%% to %.1f%%%%", p[0] * 100, p[1] * 100);
                    break;

                case kCentralInterval:
                default:
                    p[0] = 0.5 - intervals[i] / 2;
                    p[1] = 0.5 + intervals[i] / 2;
                    legend_text = "central %.1f%% interval";
                    break;
            }
            double q[2];
            if (fHistogram->GetQuantiles(2, q, p) == 2)
                hist_band = GetSubHistogram(q[0], q[1]);
        }

        if (!hist_band)
            continue;

        // set style of band histogram
        hist_band->SetFillStyle(GetBandFillStyle());
        hist_band->SetFillColor(GetBandColor(i));
        hist_band->SetLineColor(GetBandColor(i));
        hist_band->SetLineWidth(0);
        hist_band->SetLineStyle(0);
        hist_band->SetTitle(Form(legend_text.data(), hist_band->Integral("width") * 100));

        // draw band
        hist_band->Draw(options.data());
        fROOTObjects.push_back(hist_band);

    }
    for (int i = fROOTObjects.size() - 1; i >= i0; --i)
        AddBandLegendEntry(fROOTObjects[i], "", "F");

    // redraw histogram
    GetHistogram()->Draw(options.data());
}

// ---------------------------------------------------------
void BCH1D::DrawMarkers()
{
    DrawQuantiles(fNQuantiles);
    BCHistogramBase::DrawMarkers();
    DrawMedian();
}


// ---------------------------------------------------------
void BCH1D::DrawQuantiles(const unsigned n)
{
    if (n <= 1)
        return;

    // calculate quantile values (q)
    std::vector<double> q(n - 1, 0);
    std::vector<double> p(n - 1, 0);
    for (unsigned i = 1; i < n; ++i)
        p[i - 1] = i * 1. / n;
    if (GetHistogram()->GetQuantiles(n - 1, &q[0], &p[0]) != (int)n - 1)
        return;

    TLine* quantile_line = new TLine();
    fROOTObjects.push_back(quantile_line);
    quantile_line->SetLineStyle(2);
    quantile_line->SetLineColor(GetQuantileLineColor());

    double ymin = gPad->GetUymin();
    double ymax = gPad->GetUymax();
    if (gPad->GetLogy()) {
        ymin = pow(10, ymin);
        ymax = pow(10, ymax);
    }

    // draw quantile lines
    for (unsigned i = 0; i < n - 1; ++i)
        quantile_line->DrawLine(q[i], ymin, q[i], GetHistogram()->GetBinContent(GetHistogram()->FindFixBin(q[i])));

    std::string quantile_text;
    switch (n) {
        case 2:
            quantile_text = "median";
            break;
        case 3:
            quantile_text = "terciles";
            break;
        case 4:
            quantile_text = "quartiles";
            break;
        case 5:
            quantile_text = "quintiles";
            break;
        case 6:
            quantile_text = "sextiles";
            break;
        case 7:
            quantile_text = "septiles";
            break;
        case 8:
            quantile_text = "octiles";
            break;
        case 10:
            quantile_text = "deciles";
            break;
        case 12:
            quantile_text = "duodeciles";
            break;
        case 20:
            quantile_text = "vigintiles";
            break;
        case 100:
            quantile_text = "percentiles";
            break;
        default:
            quantile_text = Form("%d-quantiles", n);
            break;
    }
    AddLegendEntry(quantile_line, quantile_text.data(), "L");
}

// ---------------------------------------------------------
void BCH1D::DrawMedian()
{
    if ( !fDrawMedian || fNQuantiles == 2)
        return;

    double ymin = gPad->GetUymin();
    double ymax = gPad->GetUymax();
    double ymid = 0.5 * (ymin + ymax);
    if (gPad->GetLogy()) {
        ymin = pow(10, ymin);
        ymax = pow(10, ymax);
        ymid = pow(10, ymid);
    }

    TMarker* marker_median = new TMarker(GetMedian(), ymid * (fLogy ? pow(ymax / ymin, -0.1) : 0.8), 21);
    fROOTObjects.push_back(marker_median);
    marker_median->SetMarkerColor(GetMarkerColor());
    marker_median->SetMarkerSize(fMarkerScale * gPad->GetWNDC());
    marker_median->Draw();

    TLegendEntry* le = 0;
    double q[2], p[2] = {0.1587, 0.8413};

    if (fDrawCentral68 && GetHistogram()->GetQuantiles(2, q, p) == 2) {
        TArrow* arrow_ci = new TArrow(q[0], ymid * (fLogy ? pow(ymax / ymin, -0.1) : 0.8),
                                      q[1], ymid * (fLogy ? pow(ymax / ymin, -0.1) : 0.8),
                                      0.02 * gPad->GetWNDC(), "<|>");
        fROOTObjects.push_back(arrow_ci);
        arrow_ci->SetLineColor(marker_median->GetMarkerColor());
        arrow_ci->SetFillColor(marker_median->GetMarkerColor());
        arrow_ci->Draw();
        le = AddLegendEntry(arrow_ci, "median and central 68% interval", "PL");
        le->SetLineColor(arrow_ci->GetLineColor());
    } else
        le = AddLegendEntry(marker_median, "median", "P");
    le->SetMarkerStyle(marker_median->GetMarkerStyle());
    le->SetMarkerSize(marker_median->GetMarkerSize());
    le->SetMarkerColor(marker_median->GetMarkerColor());
}

// ---------------------------------------------------------
TH1* BCH1D::GetSubHistogram(double min, double max, const std::string& name, bool preserve_range)
{
    if (min == max || !GetHistogram())
        return 0;
    if (min > max)
        return GetSubHistogram(max, min, name);

    double xmin = GetHistogram()->GetXaxis()->GetXmin();
    double xmax = GetHistogram()->GetXaxis()->GetXmax();

    if (max < xmin || min > xmax)
        return 0;

    std::string newName(name);
    if (name.empty())
        newName = std::string(GetHistogram()->GetName()) + "_subhist";

    if ( min <= xmin && max >= xmax ) {
        return BCAux::OwnClone(GetHistogram(), name);
    }
    min = std::max<double>(min, xmin);
    max = std::min<double>(max, xmax);

    int imin = (min > xmin) ? GetHistogram()->FindFixBin(min) : 1;
    int imax = (max < xmax) ? GetHistogram()->FindFixBin(max) : GetHistogram()->GetNbinsX();

    // create binning for new histogram
    std::vector<double> bins(GetHistogram()->GetNbinsX() + 2, 0);
    bins[0] = (preserve_range) ? xmin : min;
    int i0 = (preserve_range) ? 2 : imin + 1;
    int i1 = (preserve_range) ? GetHistogram()->GetNbinsX() : imax;
    unsigned n = 1;
    for (int i = i0; i <= i1; ++i) {
        bins[n++] = GetHistogram()->GetXaxis()->GetBinLowEdge(i);
        if (min > GetHistogram()->GetXaxis()->GetBinLowEdge(i) && min < GetHistogram()->GetXaxis()->GetBinUpEdge(i))
            bins[n++] = min;
        if (max > GetHistogram()->GetXaxis()->GetBinLowEdge(i) && max < GetHistogram()->GetXaxis()->GetBinUpEdge(i))
            bins[n++] = max;
    }
    if (preserve_range || max == GetHistogram()->GetXaxis()->GetBinUpEdge(i1))
        bins[n++] = GetHistogram()->GetXaxis()->GetBinUpEdge(i1);

    // now define the new histogram
    TH1* h0;
    {
        BCAux::RootSideEffectGuard g;
        h0 = new TH1D(newName.data(), Form("%s;%s;%s", GetHistogram()->GetTitle(), GetHistogram()->GetXaxis()->GetTitle(), GetHistogram()->GetYaxis()->GetTitle()), n - 1, &bins[0]);
    }
    imin = h0->FindFixBin(min);
    imax = h0->FindFixBin(max);
    for (int i = imin; i <= imax; ++i)
        h0->SetBinContent(i, GetHistogram()->GetBinContent(GetHistogram()->FindFixBin(h0->GetBinCenter(i))));
    return h0;
}

// ---------------------------------------------------------
void BCH1D::PrintSummary(const std::string& prefix, unsigned prec, std::vector<double> intervals)
{
    if (!GetHistogram())
        return;

    double p[7] = {5e-2, 10e-2, 16e-2, 50e-2, 84e-2, 90e-2, 95e-2};
    double q[7];
    GetHistogram()->GetQuantiles(7, q, p);

    BCLog::OutSummary(prefix + Form("Mean +- sqrt(Variance):         %+.*g +- %.*g", prec, GetHistogram()->GetMean(), prec, GetHistogram()->GetRMS()));
    BCLog::OutSummary(prefix + Form("Median +- central 68%% interval: %+.*g + %.*g - %.*g", prec, q[3], prec, q[4] - q[3], prec, q[3] - q[2]));
    BCLog::OutSummary(prefix + Form("(Marginalized) mode:            %+.*g",  prec, GetLocalMode()));
    BCLog::OutSummary(prefix + Form("%2.0f%% quantile:                   %+.*g", 100 * p[0], prec, q[0]));
    BCLog::OutSummary(prefix + Form("%2.0f%% quantile:                   %+.*g", 100 * p[1], prec, q[1]));
    BCLog::OutSummary(prefix + Form("%2.0f%% quantile:                   %+.*g", 100 * p[2], prec, q[2]));
    BCLog::OutSummary(prefix + Form("%2.0f%% quantile:                   %+.*g", 100 * p[4], prec, q[4]));
    BCLog::OutSummary(prefix + Form("%2.0f%% quantile:                   %+.*g", 100 * p[5], prec, q[5]));
    BCLog::OutSummary(prefix + Form("%2.0f%% quantile:                   %+.*g", 100 * p[6], prec, q[6]));

    std::vector<BCH1D::BCH1DSmallestInterval> v = GetSmallestIntervals(intervals);
    for (unsigned i = 0; i < v.size(); ++i)
        v[i].PrintSummary(prefix, prec);
}

// ---------------------------------------------------------
void BCH1D::BCH1DSmallestInterval::PrintSummary(const std::string& prefix, unsigned prec) const
{
    BCLog::OutSummary(prefix + Form("Smallest interval%s containing %.1f%% and local mode%s:", (intervals.size() > 1 ? "s" : ""), 100. * total_mass, (intervals.size() > 1 ? "s" : "")));
    for (unsigned i = 0; i < intervals.size(); ++i)
        intervals[i].PrintSummary(prefix, prec);
}

// ---------------------------------------------------------
void BCH1D::BCH1DInterval::PrintSummary(const std::string& prefix, unsigned prec) const
{
    BCLog::OutSummary(prefix + Form("(%.*g, %.*g) (local mode at %.*g with rel. height %.*g; rel. area %.*g)",
                                    prec, xmin,
                                    prec, xmax,
                                    prec, mode,
                                    prec, relative_height,
                                    prec, relative_mass));
}

// ---------------------------------------------------------
std::vector<BCH1D::BCH1DSmallestInterval> BCH1D::GetSmallestIntervals(std::vector<double> masses)
{
    std::vector<std::pair<double, double> > bounds = GetSmallestIntervalBounds(masses);

    std::vector<BCH1D::BCH1DSmallestInterval> result;

    for (unsigned i = 0; i < bounds.size(); ++i) {
        BCH1D::BCH1DSmallestInterval smallest_interval;
        smallest_interval.total_mass = 0;
        for (int b = 1; b <= GetHistogram()->GetNbinsX(); ++b)
            if (GetHistogram()->GetBinContent(b) >= bounds[i].first) {
                BCH1D::BCH1DInterval interval;
                interval.xmin = GetHistogram()->GetXaxis()->GetBinLowEdge(b);
                interval.xmax = GetHistogram()->GetXaxis()->GetBinUpEdge(b);
                interval.relative_height = GetHistogram()->GetBinContent(b);
                interval.mode = GetHistogram()->GetBinCenter(b);
                interval.relative_mass += GetHistogram()->Integral(b, b, "width");
                while (b < GetHistogram()->GetNbinsX() and GetHistogram()->GetBinContent(b + 1) >= bounds[i].first) {
                    interval.xmax = GetHistogram()->GetXaxis()->GetBinUpEdge(++b);
                    interval.relative_mass += GetHistogram()->Integral(b, b, "width");
                    if (GetHistogram()->GetBinContent(b) > interval.relative_height) {
                        interval.relative_height = GetHistogram()->GetBinContent(b);
                        interval.mode = GetHistogram()->GetBinCenter(b);
                    }
                }
                smallest_interval.intervals.push_back(interval);
                if (smallest_interval.total_mass == 0) {
                    smallest_interval.mode = interval.mode;
                    smallest_interval.max_val = interval.relative_height;
                } else if (interval.relative_height > smallest_interval.max_val) {
                    smallest_interval.mode = interval.mode;
                    smallest_interval.max_val = interval.relative_height;
                }
                smallest_interval.total_mass += interval.relative_mass;
            }
        for (unsigned j = 0; j < smallest_interval.intervals.size(); ++j) {
            smallest_interval.intervals[j].relative_mass /= smallest_interval.total_mass;
            smallest_interval.intervals[j].relative_height /= smallest_interval.max_val;
        }
        result.push_back(smallest_interval);
    }
    return result;
}

// ---------------------------------------------------------
BCH1D::BCH1DInterval::BCH1DInterval():
    xmin(std::numeric_limits<double>::quiet_NaN()),
    xmax(std::numeric_limits<double>::quiet_NaN()),
    mode(std::numeric_limits<double>::quiet_NaN()),
    relative_height(0),
    relative_mass(0)
{
}

// ---------------------------------------------------------
BCH1D::BCH1DSmallestInterval::BCH1DSmallestInterval():
    intervals(),
    total_mass(0),
    mode(0),
    max_val(0)
{
}
