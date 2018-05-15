/*
 * Copyright (C) 2007-2013, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

// ---------------------------------------------------------

#include "BCHistogramBase.h"

#include "BCAux.h"
#include "BCLog.h"
#include "config.h"

#include <TArrow.h>
#include <TAxis.h>
#include <TH1.h>
#include <TH2D.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TMarker.h>
#include <TPad.h>
#include <TROOT.h>
#include <TString.h>

#include <algorithm>
#include <math.h>

// ---------------------------------------------------------
BCHistogramBase::BCHistogramBase(const TH1* const hist, int dimension)
    : fHistogram(0),
      fNLegendColumns(2),
      fBandOvercoverage(true),
      fBandFillStyle(1001),
      fLineColor(kBlack),
      fLineStyle(1),
      fLineWidth(1),
      fMarkerColor(kBlack),
      fMarkerScale(1.6),
      fLogx(false),
      fLogy(false),
      fLogz(false),
      fGridx(false),
      fGridy(false),
      fNBands(3),
      fNSmooth(0),
      fDrawGlobalMode(true),
      fDrawGlobalModeArrows(true),
      fGlobalModeMarkerStyle(24),
      fDrawLocalMode(false),
      fDrawLocalModeArrows(true),
      fLocalModeMarkerStyle(25),
      fDrawMean(true),
      fMeanMarkerStyle(20),
      fDrawStandardDeviation(true),
      fDrawLegend(true),
      fDrawStats(false),
      fDimension(dimension),
      fROOToptions("HIST")
{
    SetHistogram(hist);
    SetColorScheme(kGreenYellowRed);

    fLegend.SetNColumns(fNLegendColumns);
    fLegend.SetBorderSize(0);
    fLegend.SetFillColor(kWhite);
    fLegend.SetTextAlign(12);
    fLegend.SetTextSize(0.03);

    fIntervals = DefaultIntervals();
}

// ---------------------------------------------------------
BCHistogramBase::BCHistogramBase(const BCHistogramBase& other)
    : fHistogram(0)
    , fLegend(other.fLegend)
    , fGlobalMode(other.fGlobalMode)
    , fDimension(other.fDimension)
{
    SetHistogram(other.fHistogram);
    fLegend.Clear();

    CopyOptions(other);
}

// ---------------------------------------------------------
BCHistogramBase::~BCHistogramBase()
{
    delete fHistogram;
    for (unsigned i = 0; i < fROOTObjects.size(); ++i)
        delete fROOTObjects[i];
}

// ---------------------------------------------------------
void BCHistogramBase::CopyOptions(const BCHistogramBase& other)
{
    fNLegendColumns = other.fNLegendColumns;
    fBandFillStyle = other.fBandFillStyle;
    fLineColor = other.fLineColor;
    fLineStyle = other.fLineStyle;
    fLineWidth = other.fLineWidth;
    fMarkerColor = other.fMarkerColor;
    fMarkerScale = other.fMarkerScale;
    fLogx = other.fLogx;
    fLogy = other.fLogy;
    fLogz = other.fLogz;
    fGridx = other.fGridx;
    fGridy = other.fGridy;
    fNBands = other.fNBands;
    fNSmooth = other.fNSmooth;
    fDrawGlobalMode = other.fDrawGlobalMode;
    fDrawGlobalModeArrows = other.fDrawGlobalModeArrows;
    fGlobalModeMarkerStyle = other.fGlobalModeMarkerStyle;
    fDrawLocalMode = other.fDrawLocalMode;
    fDrawLocalModeArrows = other.fDrawLocalModeArrows;
    fLocalModeMarkerStyle = other.fLocalModeMarkerStyle;
    fDrawMean = other.fDrawMean;
    fMeanMarkerStyle = other.fMeanMarkerStyle;
    fDrawStandardDeviation = other.fDrawStandardDeviation;
    fDrawLegend = other.fDrawLegend;
    fDrawStats = other.fDrawStats;
    fBandColors = other.fBandColors;
    fBandOvercoverage = other.fBandOvercoverage;
    fIntervals = other.fIntervals;
    fROOToptions = other.fROOToptions;
}

// ---------------------------------------------------------
BCHistogramBase& BCHistogramBase::operator=(BCHistogramBase other)
{
    // swap into this
    swap(*this, other);

    return *this;
}

// ---------------------------------------------------------
void swap(BCHistogramBase& first, BCHistogramBase& second)
{
    // swap histogram pointers
    std::swap(first.fHistogram,             second.fHistogram);

    // swap modes
    std::swap(first.fLocalMode,             second.fLocalMode);
    std::swap(first.fGlobalMode,            second.fGlobalMode);

    // swap legends
    TLegend temp(first.fLegend);
    second.fLegend.Copy(first.fLegend);
    temp.Copy(second.fLegend);
    // std::swap(first.fLegend,                second.fLegend);

    // swap dimension
    std::swap(first.fDimension,             second.fDimension);

    // swap ROOT object pointer vectors
    std::swap(first.fROOTObjects,           second.fROOTObjects);

    // swap options
    std::swap(first.fNLegendColumns,        second.fNLegendColumns);
    std::swap(first.fBandFillStyle,         second.fBandFillStyle);
    std::swap(first.fLineColor,             second.fLineColor);
    std::swap(first.fLineStyle,             second.fLineStyle);
    std::swap(first.fLineWidth,             second.fLineWidth);
    std::swap(first.fMarkerColor,           second.fMarkerColor);
    std::swap(first.fMarkerScale,           second.fMarkerScale);
    std::swap(first.fLogx,                  second.fLogx);
    std::swap(first.fLogy,                  second.fLogy);
    std::swap(first.fLogz,                  second.fLogz);
    std::swap(first.fGridx,                 second.fGridx);
    std::swap(first.fGridy,                 second.fGridy);
    std::swap(first.fNBands,                second.fNBands);
    std::swap(first.fNSmooth,               second.fNSmooth);
    std::swap(first.fDrawGlobalMode,        second.fDrawGlobalMode);
    std::swap(first.fDrawGlobalModeArrows,  second.fDrawGlobalModeArrows);
    std::swap(first.fGlobalModeMarkerStyle, second.fGlobalModeMarkerStyle);
    std::swap(first.fDrawLocalMode,         second.fDrawLocalMode);
    std::swap(first.fDrawLocalModeArrows,   second.fDrawLocalModeArrows);
    std::swap(first.fLocalModeMarkerStyle,  second.fLocalModeMarkerStyle);
    std::swap(first.fDrawMean,              second.fDrawMean);
    std::swap(first.fMeanMarkerStyle,       second.fMeanMarkerStyle);
    std::swap(first.fDrawStandardDeviation, second.fDrawStandardDeviation);
    std::swap(first.fDrawLegend,            second.fDrawLegend);
    std::swap(first.fDrawStats,             second.fDrawStats);
    std::swap(first.fBandColors,            second.fBandColors);
    std::swap(first.fBandOvercoverage,      second.fBandOvercoverage);
    std::swap(first.fIntervals,             second.fIntervals);
    std::swap(first.fROOToptions,           second.fROOToptions);
}

// ---------------------------------------------------------
void BCHistogramBase::SetHistogram(const TH1* const hist)
{
    delete fHistogram;

    if (!hist || (fDimension >= 0 && hist->GetDimension() != fDimension)) {
        fHistogram = 0;
        fLocalMode.clear();
        return;
    }

    fHistogram = BCAux::OwnClone(hist, Form("%s_bch", hist->GetName()));
    fHistogram->SetStats(false);
    fHistogram->SetDirectory(0);
    fDimension = fHistogram->GetDimension();

    // normalize; TODO: replace with division of each bin by width/area for arbitrary binning
    double integral = GetHistogram()->Integral("width");
    if (integral != 0)
        GetHistogram()->Scale(1. / integral);

    // Get local mode
    int b = GetHistogram()->GetMaximumBin();
    int bx, by, bz;
    GetHistogram()->GetBinXYZ(b, bx, by, bz);
    fLocalMode.assign(1, GetHistogram()->GetXaxis()->GetBinCenter(bx));
    if (by > 0)
        fLocalMode.push_back(GetHistogram()->GetYaxis()->GetBinCenter(by));
    if (bz > 0)
        fLocalMode.push_back(GetHistogram()->GetZaxis()->GetBinCenter(bz));
}

// ---------------------------------------------------------
bool BCHistogramBase::Valid() const
{
    return fHistogram and (fHistogram->Integral() != 0);
}

// ---------------------------------------------------------
void BCHistogramBase::SetColorScheme(BCHColorScheme scheme)
{
    fBandColors.clear();

    switch (scheme) {

        case kBlackWhite:
            AddBandColor(12);						// dark
            AddBandColor(14);
            AddBandColor(16);
            AddBandColor(17);						// to
            AddBandColor(18);
            AddBandColor(19);						// light
            AddBandColor(10);						// white
            SetMarkerColor(kBlack);
            SetLineColor(kBlack);
            break;

        case kBlueOrange:
            AddBandColor(kBlue);
            AddBandColor(kBlue - 3);
            AddBandColor(kBlue - 1);
            AddBandColor(kBlue - 6);
            AddBandColor(kBlue - 8);
            AddBandColor(kBlue - 9);
            AddBandColor(kBlue - 10);
            SetMarkerColor(kOrange);
            SetLineColor(kBlack);
            break;

        case kRedGreen:
            AddBandColor(kRed);
            AddBandColor(kRed - 3);
            AddBandColor(kRed - 1);
            AddBandColor(kRed - 6);
            AddBandColor(kRed - 8);
            AddBandColor(kRed - 9);
            AddBandColor(kRed - 10);
            SetMarkerColor(kGreen);
            SetLineColor(kBlack);
            break;

        case kGreenYellowRed:
        default:
            AddBandColor(kGreen);
            AddBandColor(kYellow);
            AddBandColor(kRed);
            AddBandColor(kRed - 3);
            AddBandColor(kRed - 1);
            AddBandColor(kRed - 6);
            SetMarkerColor(kBlack);
            SetLineColor(kBlack);
            break;

    }

}

// ---------------------------------------------------------
void BCHistogramBase::Smooth(int n)
{
    if (n < 0)
        n = fNSmooth;
    if (n <= 0)
        return;
    GetHistogram()->Scale(GetHistogram()->Integral("width"));
    if (GetHistogram()->GetDimension() == 1) {
        if (GetHistogram()->GetNbinsX() >= 3)
            GetHistogram()->Smooth(n);
    } else // ROOT currently only allows single instances of smoothing for TH2
        for (int i = 0; i < n; ++i)
            GetHistogram()->Smooth(1);
    double integral = GetHistogram()->Integral("width");
    if ( integral != 0 )
        GetHistogram()->Scale(1. / integral);
}

// ---------------------------------------------------------
void BCHistogramBase::CheckIntervals(std::vector<double>& intervals, int sort)
{
    // remove out-of-bounds entries
    for (int i = intervals.size() - 1; i >= 0; --i)
        if (intervals[i] < 0 || intervals[i] > 1) {
            BCLog::OutWarning(Form("BCHistogramBase::CheckIntervals : interval out of bounds, removing %f", intervals[i]));
            intervals.erase(intervals.begin() + i);
        }

    // sort
    if (sort > 0)
        std::stable_sort(intervals.begin(), intervals.end(), std::less<double>());
    if (sort < 0)
        std::stable_sort(intervals.begin(), intervals.end(), std::greater<double>());

    // warning: unique will not reliably remove all duplicates unless the vector has first been sorted.
    unsigned old_size = intervals.size();
    std::vector<double>::iterator last = std::unique(intervals.begin(), intervals.end());
    intervals.erase(last, intervals.end());
    if (intervals.size() < old_size)
        BCLog::OutWarning(Form("BCHistogramBase::CheckIntervals : %lu duplicate interval values were removed.", old_size - intervals.size()));
}

// ---------------------------------------------------------
std::vector<double> BCHistogramBase::DefaultIntervals(int nbands)
{
    if (nbands < 0)
        nbands = fNBands;
    std::vector<double> intervals;
    if (nbands > 0) intervals.push_back(0.682689492137);
    if (nbands > 1) intervals.push_back(0.954499736104);
    if (nbands > 2) intervals.push_back(0.997300203937);
    if (nbands > 3) intervals.push_back(0.999936657516);
    if (nbands > 4) intervals.push_back(0.999999426697);
    if (nbands > 5) intervals.push_back(0.999999998027);
    if (nbands > 6) intervals.push_back(0.999999999997);
    if (nbands > 7) intervals.push_back(1);
    return intervals;
}

// ---------------------------------------------------------
void BCHistogramBase::GetNonzeroBinDensityMassVector(std::vector<std::pair<double, double> >& bin_dens_mass, int sort)
{
    if (!GetHistogram())
        return;

    // calculate number of bins
    unsigned nbins = GetHistogram()->GetBin(GetHistogram()->GetNbinsX(), GetHistogram()->GetNbinsY(), GetHistogram()->GetNbinsZ());

    // reserve space
    bin_dens_mass.reserve(nbins);

    // fill bin_dens_mass with pairs of (prob. density, prob. mass) for all non-empty bins
    for (unsigned i = 1; i <= nbins; ++i)
        if (!GetHistogram()->IsBinUnderflow(i) && !GetHistogram()->IsBinOverflow(i) && GetHistogram()->GetBinContent(i) > 0) {
            // if 1D, by = bz = -1; if 2D, bz = -1
            int bx, by, bz;
            GetHistogram()->GetBinXYZ(i, bx, by, bz);
            // bin val
            double v = GetHistogram()->GetBinContent(i);
            // bin width / area / volume
            double w = (bx > 0 ? GetHistogram()->GetXaxis()->GetBinWidth(bx) : 1) * (by > 0 ? GetHistogram()->GetYaxis()->GetBinWidth(by) : 1) * (bz > 0 ? GetHistogram()->GetZaxis()->GetBinWidth(bz) : 1);
            bin_dens_mass.push_back(std::make_pair(v, v * w));
        }
    // sort
    if (sort < 0)
        std::stable_sort(bin_dens_mass.begin(), bin_dens_mass.end(), std::greater<std::pair<double, double> >());
    if (sort > 0)
        std::stable_sort(bin_dens_mass.begin(), bin_dens_mass.end(), std::less<std::pair<double, double> >());
}

// ---------------------------------------------------------
std::vector<std::pair<double, double> > BCHistogramBase::GetSmallestIntervalBounds(std::vector<double> masses, bool overcoverage)
{
    std::vector<std::pair<double, double> > levels;

    // if no masses asked for or no histogram set, return empty vector
    if (masses.empty() || !GetHistogram())
        return levels;

    // check and sort masses to decreasing order
    CheckIntervals(masses, -1);

    // reserce space in the levels vector
    levels.reserve(masses.size());

    std::vector<std::pair<double, double> > bin_dens_mass;
    GetNonzeroBinDensityMassVector(bin_dens_mass, -1); // sorted in decreasing order

    // consolidate bins by density
    std::vector<std::pair<double, double> > dens_mass;
    dens_mass.reserve(bin_dens_mass.size());
    for (unsigned i = 0; i < bin_dens_mass.size(); ++i) {
        dens_mass.push_back(bin_dens_mass[i]);
        while (i < bin_dens_mass.size() - 1 and bin_dens_mass[i + 1].first == dens_mass.back().first)
            dens_mass.back().second += bin_dens_mass[++i].second;
    }

    if (overcoverage) {
        double prob_sum = 0;
        for (unsigned i = 0; i < dens_mass.size() and !masses.empty() and prob_sum <= 1; ++i) {
            prob_sum += dens_mass[i].second;
            while (!masses.empty() and prob_sum >= masses.back()) {
                levels.push_back(std::make_pair(dens_mass[i].first, prob_sum));
                masses.pop_back();
            }
        }
    }	else { // undercoverage
        // remove mass levels that cannot be undercovered
        double prob_sum = dens_mass.front().second;
        while (!masses.empty() and masses.back() < prob_sum)
            masses.pop_back();
        for (unsigned i = 1; i < dens_mass.size() and !masses.empty() and prob_sum <= 1; ++i) {
            prob_sum += dens_mass[i].second;
            while (!masses.empty() and prob_sum >= masses.back()) {
                levels.push_back(std::make_pair(dens_mass[i - 1].first, prob_sum - dens_mass[i].second));
                masses.pop_back();
            }
        }
    }

    std::vector<std::pair<double, double> >::iterator last = std::unique(levels.begin(), levels.end());
    levels.erase(last, levels.end());
    return levels;
}

// ---------------------------------------------------------
std::vector<double> BCHistogramBase::GetSmallestIntervalSize(std::vector<double> masses, bool overcoverage)
{
    // vector of sizes
    std::vector<double> S;

    // remove out-of-bounds entries
    for (int i = masses.size() - 1; i >= 0; --i)
        if (masses[i] < 0 || masses[i] > 1) {
            BCLog::OutWarning(Form("BCHistogramBase::GetSmallestIntervalSize : mass out of bounds, removing %f", masses[i]));
            masses.erase(masses.begin() + i);
        }
    if (masses.empty())
        return S;

    std::vector<std::pair<double, double> > bin_dens_mass;
    GetNonzeroBinDensityMassVector(bin_dens_mass, -1); // sorted in decreasing order

    if (bin_dens_mass.empty())
        return S;

    S.assign(masses.size(), 0);
    double mass = 0;
    for (std::vector<std::pair<double, double> >::const_iterator bin = bin_dens_mass.begin(); bin != bin_dens_mass.end(); ++bin) {
        for (unsigned i = 0; i < masses.size(); ++i)
            // if more mass is needed, add area
            if (mass + (overcoverage ? 0 : bin->second) < masses[i])
                S[i] += bin->second / bin->first;
        mass += bin->second;
    }
    return S;
}

// ---------------------------------------------------------
double BCHistogramBase::GetSmallestIntervalSize(double mass, bool overcoverage)
{
    std::vector<double> s = GetSmallestIntervalSize(std::vector<double>(1, mass), overcoverage);
    if (s.empty())
        return 0;
    return s[0];
}

// ---------------------------------------------------------
void BCHistogramBase::Draw()
{
    if (!GetHistogram())
        return;

    GetHistogram()->SetStats(fDrawStats);
    fLegend.SetNColumns(fNLegendColumns);

    Smooth(fNSmooth);

    std::string options = fROOToptions;
    std::transform(options.begin(), options.end(), options.begin(), ::tolower);
    // if option "same" is not specified, draw axes and add "same" to options
    if (options.find("same") == std::string::npos) {
        gPad->SetLogx(fLogx);
        gPad->SetLogy(fLogy);
        gPad->SetLogz(fLogz);
        gPad->SetGridx(fGridx);
        gPad->SetGridy(fGridy);
        if (GetHistogram()->GetDimension() == 1) {
            // necessary because ROOT will otherwise draw space below zero for the errors---which we don't have
            double ymin = GetHistogram()->GetMinimum();
            double ymax = GetHistogram()->GetMaximum();
            GetHistogram()->GetYaxis()->SetRangeUser(ymin, (gPad->GetLogy() ? 2 * ymax : ymax + 0.1 * (ymax - ymin)));
        }
        GetHistogram()->Draw("axis");
        options += "same";
    }

    GetHistogram()->SetLineColor(GetLineColor());
    GetHistogram()->SetLineStyle(GetLineStyle());
    GetHistogram()->SetLineWidth(GetLineWidth());
    DrawBands(options);

    DrawMarkers();
    if (fDrawLegend)
        DrawLegend();

    gPad->RedrawAxis();
    gPad->Update();
}

// ---------------------------------------------------------
void BCHistogramBase::DrawMarkers()
{
    DrawGlobalMode();
    DrawLocalMode();
    DrawMean();
}

// ---------------------------------------------------------
void BCHistogramBase::DrawGlobalMode()
{
    gPad->Update();
    double ymin = gPad->GetUymin();
    double ymax = gPad->GetUymax();
    double y = ymin + 0.3 * (ymax - ymin);
    if (gPad->GetLogy()) {
        ymin = pow(10, ymin);
        ymax = pow(10, ymax);
        y = ymin * pow(ymax / ymin, 0.3);
    }
    if (GetHistogram()->GetDimension() > 1 && fGlobalMode.size() > 1)
        y = fGlobalMode[1];

    if (fDrawGlobalMode && !fGlobalMode.empty()) {
        TMarker* marker_mode = new TMarker(fGlobalMode[0], y, fGlobalModeMarkerStyle);
        fROOTObjects.push_back(marker_mode);
        marker_mode->SetMarkerColor(GetMarkerColor());
        marker_mode->SetMarkerSize(fMarkerScale * gPad->GetWNDC());
        marker_mode->Draw();

        TLegendEntry* le = AddLegendEntry(marker_mode, "global mode", "P");
        le->SetMarkerStyle(marker_mode->GetMarkerStyle());
        le->SetMarkerSize(marker_mode->GetMarkerSize());
        le->SetMarkerColor(marker_mode->GetMarkerColor());

        if (fDrawGlobalModeArrows) {
            TArrow* arrow_mode = new TArrow(marker_mode->GetX(), (gPad->GetLogy() ? marker_mode->GetY()*pow(ymax / ymin, -1.5e-2) : marker_mode->GetY() + (ymax - ymin) * -1.5e-2),
                                            marker_mode->GetX(), (gPad->GetLogy() ? ymin * pow(ymax / ymin, 4e-2) : ymin + (ymax - ymin) * 4e-2),
                                            2e-2 * gPad->GetWNDC(), "|>");
            fROOTObjects.push_back(arrow_mode);
            arrow_mode->SetLineColor(marker_mode->GetMarkerColor());
            arrow_mode->SetFillColor(marker_mode->GetMarkerColor());
            arrow_mode->Draw();

            if (GetHistogram()->GetDimension() > 1 && fGlobalMode.size() > 1) {
                double xmin = gPad->GetUxmin();
                double xmax = gPad->GetUxmax();
                if (gPad->GetLogx()) {
                    ymin = pow(10, xmin);
                    ymax = pow(10, xmax);
                }
                TArrow* arrow_mode2 = new TArrow((gPad->GetLogx() ? marker_mode->GetX()*pow(xmax / xmin, -1.5e-2) : marker_mode->GetX() + (xmax - xmin) * -1.5e-2), marker_mode->GetY(),
                                                 (gPad->GetLogx() ? xmin * pow(xmax / xmin, 4e-2) : xmin + (xmax - xmin) * 4e-2), marker_mode->GetY(),
                                                 2e-2 * gPad->GetWNDC(), "|>");
                fROOTObjects.push_back(arrow_mode2);
                arrow_mode2->SetLineColor(marker_mode->GetMarkerColor());
                arrow_mode2->SetFillColor(marker_mode->GetMarkerColor());
                arrow_mode2->Draw();
            }
        }

    }
}

// ---------------------------------------------------------
void BCHistogramBase::DrawLocalMode()
{
    gPad->Update();
    double ymin = gPad->GetUymin();
    double ymax = gPad->GetUymax();
    double y = ymin + 0.25 * (ymax - ymin);
    if (gPad->GetLogy()) {
        ymin = pow(10, ymin);
        ymax = pow(10, ymax);
        y = ymin * pow(ymax / ymin, 0.25);
    }
    if (GetHistogram()->GetDimension() > 1 && fLocalMode.size() > 1)
        y = fLocalMode[1];

    if (fDrawLocalMode && !fLocalMode.empty()) {
        TMarker* marker_mode = new TMarker(fLocalMode[0], y, fLocalModeMarkerStyle);
        fROOTObjects.push_back(marker_mode);
        marker_mode->SetMarkerColor(GetMarkerColor());
        marker_mode->SetMarkerSize(fMarkerScale * gPad->GetWNDC());
        marker_mode->Draw();

        TLegendEntry* le = AddLegendEntry(marker_mode, "local mode", "P");
        le->SetMarkerStyle(marker_mode->GetMarkerStyle());
        le->SetMarkerSize(marker_mode->GetMarkerSize());
        le->SetMarkerColor(marker_mode->GetMarkerColor());

        if (fDrawLocalModeArrows) {
            TArrow* arrow_mode = new TArrow(marker_mode->GetX(), (gPad->GetLogy() ? marker_mode->GetY()*pow(ymax / ymin, -1.5e-2) : marker_mode->GetY() + (ymax - ymin) * -1.5e-2),
                                            marker_mode->GetX(), (gPad->GetLogy() ? ymin * pow(ymax / ymin, 3e-2) : ymin + (ymax - ymin) * 3e-2),
                                            2e-2 * gPad->GetWNDC(), "|>");
            fROOTObjects.push_back(arrow_mode);
            arrow_mode->SetLineColor(marker_mode->GetMarkerColor());
            arrow_mode->SetFillColor(marker_mode->GetMarkerColor());
            arrow_mode->Draw();

            if (GetHistogram()->GetDimension() > 1 && fLocalMode.size() > 1) {
                double xmin = gPad->GetUxmin();
                double xmax = gPad->GetUxmax();
                if (gPad->GetLogx()) {
                    ymin = pow(10, xmin);
                    ymax = pow(10, xmax);
                }
                TArrow* arrow_mode2 = new TArrow((gPad->GetLogx() ? marker_mode->GetX()*pow(xmax / xmin, -1.5e-2) : marker_mode->GetX() + (xmax - xmin) * -1.5e-2), marker_mode->GetY(),
                                                 (gPad->GetLogx() ? xmin * pow(xmax / xmin, 3e-2) : xmin + (xmax - xmin) * 3e-2), marker_mode->GetY(),
                                                 2e-2 * gPad->GetWNDC(), "|>");
                fROOTObjects.push_back(arrow_mode2);
                arrow_mode2->SetLineColor(marker_mode->GetMarkerColor());
                arrow_mode2->SetFillColor(marker_mode->GetMarkerColor());
                arrow_mode2->Draw();
            }
        }

    }
}

// ---------------------------------------------------------
void BCHistogramBase::DrawMean()
{
    gPad->Update();
    double ymin = gPad->GetUymin();
    double ymax = gPad->GetUymax();
    double y = ymin + 0.4 * (ymax - ymin);
    if (gPad->GetLogy()) {
        ymin = pow(10, ymin);
        ymax = pow(10, ymax);
        y = ymin * pow(ymax / ymin, 0.4);
    }
    if (GetHistogram()->GetDimension() > 1)
        y = GetHistogram()->GetMean(2);

    if ( fDrawMean ) {
        TMarker* marker_mean = new TMarker(GetHistogram()->GetMean(1), y, fMeanMarkerStyle);
        fROOTObjects.push_back(marker_mean);
        marker_mean->SetMarkerColor(GetMarkerColor());
        marker_mean->SetMarkerSize(fMarkerScale * gPad->GetWNDC());
        marker_mean->Draw();

        // legend entry is managed separately and need not be in trash
        TLegendEntry* le = 0;
        if ( fDrawStandardDeviation ) {
            TArrow* arrow_std = new TArrow(marker_mean->GetX() - GetHistogram()->GetRMS(1), marker_mean->GetY(),
                                           marker_mean->GetX() + GetHistogram()->GetRMS(1), marker_mean->GetY(),
                                           0.02 * gPad->GetWNDC(), "<|>");
            fROOTObjects.push_back(arrow_std);
            arrow_std->SetLineColor(marker_mean->GetMarkerColor());
            arrow_std->SetFillColor(marker_mean->GetMarkerColor());
            arrow_std->Draw();
            le = AddLegendEntry(arrow_std, "mean and std. dev.", "PL");
            le->SetLineColor(arrow_std->GetLineColor());

            if (GetHistogram()->GetDimension() > 1) {
                TArrow* arrow_std2 = new TArrow(marker_mean->GetX(), marker_mean->GetY() - GetHistogram()->GetRMS(2),
                                                marker_mean->GetX(), marker_mean->GetY() + GetHistogram()->GetRMS(2),
                                                0.02 * gPad->GetWNDC(), "<|>");
                fROOTObjects.push_back(arrow_std2);
                arrow_std2->SetLineColor(marker_mean->GetMarkerColor());
                arrow_std2->SetFillColor(marker_mean->GetMarkerColor());
                arrow_std2->Draw();
            }
        } else {
            le = AddLegendEntry(marker_mean, "mean", "P");
        }
        le->SetMarkerStyle(marker_mean->GetMarkerStyle());
        le->SetMarkerSize(marker_mean->GetMarkerSize());
        le->SetMarkerColor(marker_mean->GetMarkerColor());
    }
}

// ---------------------------------------------------------
double BCHistogramBase::ResizeLegend()
{
#if ROOTVERSION >= 6000000
    fLegend.SetX1(gPad->GetLeftMargin() + (1 - gPad->GetRightMargin() - gPad->GetLeftMargin()) * 10e-2);
    fLegend.SetX2(1 - gPad->GetRightMargin());
    fLegend.SetY1(1 - gPad->GetTopMargin() - fLegend.GetTextSize()*fLegend.GetNRows());
    fLegend.SetY2(1 - gPad->GetTopMargin());
    fLegend.SetColumnSeparation(0.0);
    return fLegend.GetY1();
#else
    fLegend.SetX1NDC(gPad->GetLeftMargin() + (1 - gPad->GetRightMargin() - gPad->GetLeftMargin()) * 10e-2);
    fLegend.SetX2NDC(1 - gPad->GetRightMargin());
    fLegend.SetY1NDC(1 - gPad->GetTopMargin() - fLegend.GetTextSize()*fLegend.GetNRows());
    fLegend.SetY2NDC(1 - gPad->GetTopMargin());
    fLegend.SetColumnSeparation(0.0);
    return fLegend.GetY1NDC();
#endif
}

// ---------------------------------------------------------
TLegendEntry* BCHistogramBase::AddLegendEntry(TObject* obj, const std::string& label, const std::string& options)
{
    if (fExtraLegendEntries.empty())
        return fLegend.AddEntry(obj, label.data(), options.data());
    TLegendEntry* le = fExtraLegendEntries.front();
    le->SetObject(obj);
    le->SetLabel(label.data());
    le->SetOption(options.data());
    fExtraLegendEntries.erase(fExtraLegendEntries.begin());
    return le;
}

// ---------------------------------------------------------
TLegendEntry* BCHistogramBase::AddBandLegendEntry(TObject* obj, const std::string& label, const std::string& options)
{
    TLegendEntry* le = fLegend.AddEntry(obj, label.data(), options.data());
    for (int i = 1; i < fLegend.GetNColumns(); ++i)
        fExtraLegendEntries.push_back(fLegend.AddEntry((TObject*)0, "", ""));
    return le;
}

// ---------------------------------------------------------
void BCHistogramBase::DrawLegend()
{
    double ymin = gPad->GetUymin();
    double ymax = gPad->GetUymax();
    if (gPad->GetLogy()) {
        ymin = pow(10, ymin);
        ymax = pow(10, ymax);
    }

    fHistogram->GetYaxis()->SetRangeUser(ymin, ymax * (1.15 + fLegend.GetTextSize()*fLegend.GetNRows()) * 1.05);

    gPad->SetTopMargin(0.02);

    double y1ndc = ResizeLegend();
    fLegend.Draw();

    // rescale top margin
    gPad->SetTopMargin(1 - y1ndc + 0.01);

}
