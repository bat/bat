/*
 * Copyright (C) 2007-2014, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

// ---------------------------------------------------------

#include "BCH2D.h"
#include "BCH1D.h"

#include "BCLog.h"
#include "BCMath.h"

#include <TArrow.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TH2.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TMarker.h>
#include <TObject.h>
#include <TROOT.h>
#include <TString.h>
#include <TStyle.h>
#include <TLine.h>

#include <math.h>
#include <algorithm>

// ---------------------------------------------------------
BCH2D::BCH2D(const TH2* const h)
    : BCHistogramBase(h, 2)
    , fBandType(kSmallestInterval)
    , fDrawProfileX(false)
    , fProfileXType(kProfileMean)
    , fProfileXLineColor(kBlack)
    , fProfileXLineStyle(2)
    , fDrawProfileY(false)
    , fProfileYType(kProfileMean)
    , fProfileYLineColor(kBlack)
    , fProfileYLineStyle(2)
{
    SetDrawLocalMode(true, false);
    SetDrawGlobalMode(true, false);
}

// ---------------------------------------------------------
BCH2D::BCH2D(const BCH2D& other)
    : BCHistogramBase(other)
{
    CopyOptions(other);
}

// ---------------------------------------------------------
BCH2D::~BCH2D()
{
}

// ---------------------------------------------------------
void BCH2D::CopyOptions(const BCH2D& other)
{
    BCHistogramBase::CopyOptions(other);
    fBandType = other.fBandType;
    fDrawProfileX = other.fDrawProfileX;
    fProfileXType = other.fProfileXType;
    fProfileXLineColor = other.fProfileXLineColor;
    fProfileXLineStyle = other.fProfileXLineStyle;
    fDrawProfileY = other.fDrawProfileY;
    fProfileYType = other.fProfileYType;
    fProfileYLineColor = other.fProfileYLineColor;
    fProfileYLineStyle = other.fProfileYLineStyle;
}

// ---------------------------------------------------------
BCH2D& BCH2D::operator=(const BCH2D& other)
{
    // copy other
    BCH2D temp(other);
    // swap into this
    swap(*this, temp);
    return *this;
}

// ---------------------------------------------------------
void swap(BCH2D& first, BCH2D& second)
{
    // swap base class
    swap(static_cast<BCHistogramBase&>(first), static_cast<BCHistogramBase&>(second));
    // swap remaining members
    std::swap(first.fBandType,          second.fBandType);
    std::swap(first.fDrawProfileX,      second.fDrawProfileX);
    std::swap(first.fProfileXType,      second.fProfileXType);
    std::swap(first.fProfileXLineColor, second.fProfileXLineColor);
    std::swap(first.fProfileXLineStyle, second.fProfileXLineStyle);
    std::swap(first.fDrawProfileY,      second.fDrawProfileY);
    std::swap(first.fProfileYType,      second.fProfileYType);
    std::swap(first.fProfileYLineColor, second.fProfileYLineColor);
    std::swap(first.fProfileYLineStyle, second.fProfileYLineStyle);
}

// ---------------------------------------------------------
void BCH2D::CheckIntervals(std::vector<double>& intervals)
{
    if (fBandType != kNoBands)
        BCHistogramBase::CheckIntervals(intervals, +1);
}

// ---------------------------------------------------------
std::vector<double> BCH2D::DefaultIntervals(int nbands)
{
    std::vector<double> intervals;

    switch (fBandType) {

        case kNoBands:
            return intervals;

        case kSmallestInterval:
        default:
            return BCHistogramBase::DefaultIntervals(nbands);

    }
}

// ---------------------------------------------------------
void BCH2D::DrawBands(std::string options)
{
    if (fBandType == kNoBands) {
        GetHistogram()->Draw((options + "colz").data());
        gPad->Update();
        return;
    }

    std::vector<double> intervals = fIntervals;
    CheckIntervals(intervals);

    if (intervals.empty())
        return;

    // set contour levels
    std::vector<double> levels;
    std::vector<std::string> legend_text;
    switch (fBandType) {

        // 	// using 1D slicings
        // case kCentralIntervalOfYGivenX:
        // case kCentralIntervalOfXGivenY:
        // case kSmallestIntervalOfYGivenX:
        // case kSmallestIntervalOfXGivenY:
        // case kUpperLimitOfYGivenX:
        // case kUpperLimitOfXGivenY:
        // case kLowerLimitOfYGivenX:
        // case kLowerLimitOfXGivenY:
        // 	break;

        case kSmallestInterval:
        default:
            std::vector<std::pair<double, double> > dens_mass = GetSmallestIntervalBounds(intervals, fBandOvercoverage);
            for (unsigned i = 0; i < dens_mass.size(); ++i) {
                levels.push_back(dens_mass[dens_mass.size() - i - 1].first);
                legend_text.push_back(Form("smallest %.1f %% interval(s)", 100 * dens_mass[i].second));
            }
            // levels.push_back(GetHistogram()->GetMaximum());
            break;
    }

    // make sure enough colors have been designated
    while (levels.size() > fBandColors.size())
        fBandColors.push_back(fBandColors.back() - 1);

    // set contour colors
    std::vector<int> colors;
    for (int i = levels.size() - 1; i >= 0; --i)
        colors.push_back(fBandColors[i]);

    // set contour levels
    GetHistogram()->SetContour(levels.size(), &levels[0]);

    if (fBandFillStyle <= 0) {
        GetHistogram()->SetLineColor(GetLineColor());
        GetHistogram()->Draw((options + "cont2").data());
    } else {
        gStyle->SetPalette(colors.size(), &colors[0]);
        GetHistogram()->SetFillStyle(fBandFillStyle);
        GetHistogram()->Draw((options + "cont").data());
    }
    gPad->Update();

    // Set legend entries
    for (unsigned i = 0; i < levels.size(); ++i) {
        if (fBandFillStyle > 0) {
            TLegendEntry* le = AddBandLegendEntry((TObject*)0, legend_text[i].data(), "F");
            le->SetFillColor(colors[levels.size() - i - 1]);
            le->SetFillStyle(1001/*fBandFillStyle*/);
            le->SetLineColor(0);
            le->SetLineWidth(0);
            le->SetLineStyle(0);
        } else {
            TLegendEntry* le = AddBandLegendEntry((TObject*)0, legend_text[i].data(), "L");
            le->SetLineColor(GetLineColor());
            le->SetLineStyle(levels.size() - i);
        }
    }
}

// ---------------------------------------------------------
void BCH2D::DrawMarkers()
{
    BCHistogramBase::DrawMarkers();
    DrawProfileGraphs();
}

// ---------------------------------------------------------
TGraph* BCH2D::CalculateProfileGraph(BCH2DProfileAxis axis, BCH2DProfileType bt)
{

    unsigned n_i = (axis == kProfileY) ? GetHistogram()->GetNbinsY() : GetHistogram()->GetNbinsX();
    unsigned n_j = (axis == kProfileY) ? GetHistogram()->GetNbinsX() : GetHistogram()->GetNbinsY();

    TGraph* g = new TGraph();

    // loop over bins of axis to be profiled
    for (unsigned i = 1; i <= n_i; ++i) {

        switch (bt) {

            case kProfileMedian: {

                // calculate 50% of total probability mass in slice
                double median_prob = 0.5 * ((axis == kProfileY) ? GetHistogram()->Integral(1, n_j, i, i, "width") : GetHistogram()->Integral(i, i, 1, n_j, "width"));
                if (median_prob <= 0)
                    break;

                double prob_sum = 0;
                // loop until 50% of probability mass is found
                for (unsigned j = 1; j <= n_j; ++j) {
                    prob_sum += (axis == kProfileY) ? GetHistogram()->Integral(j, j, i, i, "width") : GetHistogram()->Integral(i, i, j, j, "width");
                    if (prob_sum > median_prob) {
                        if (axis == kProfileY)
                            g->SetPoint(g->GetN(), GetHistogram()->GetXaxis()->GetBinLowEdge(j), GetHistogram()->GetYaxis()->GetBinCenter(i));
                        else
                            g->SetPoint(g->GetN(), GetHistogram()->GetXaxis()->GetBinCenter(i), GetHistogram()->GetYaxis()->GetBinLowEdge(j));
                        break;
                    }
                }
                break;
            }

            case kProfileMode: {
                double max_val = 0;
                unsigned j_max_val = 0;
                for (unsigned j = 1; j <= n_j; ++j) {
                    double val = (axis == kProfileY) ? GetHistogram()->GetBinContent(j, i) : GetHistogram()->GetBinContent(i, j);
                    if (val > max_val) {
                        j_max_val = j;
                        max_val = val;
                    }
                }
                if (j_max_val > 0) {
                    if (axis == kProfileY)
                        g->SetPoint(g->GetN(), GetHistogram()->GetXaxis()->GetBinCenter(j_max_val), GetHistogram()->GetYaxis()->GetBinCenter(i));
                    else
                        g->SetPoint(g->GetN(), GetHistogram()->GetXaxis()->GetBinCenter(i), GetHistogram()->GetYaxis()->GetBinCenter(j_max_val));
                }
                break;
            }

            case kProfileMean:
            default: {
                // calculate total probability mass in slice
                double mass_sum = 0;
                double sum = 0 ;
                for (unsigned j = 0; j <= n_j; ++j) {
                    double mass = (axis == kProfileY) ? GetHistogram()->Integral(j, j, i, i, "width") : GetHistogram()->Integral(i, i, j, j, "width");
                    mass_sum += mass;
                    sum += mass * ((axis == kProfileY) ? GetHistogram()->GetXaxis()->GetBinCenter(j) : GetHistogram()->GetYaxis()->GetBinCenter(j));
                }
                if (mass_sum >= 0) {
                    if (axis == kProfileY)
                        g->SetPoint(g->GetN(), sum / mass_sum, GetHistogram()->GetYaxis()->GetBinCenter(i));
                    else
                        g->SetPoint(g->GetN(), GetHistogram()->GetXaxis()->GetBinCenter(i), sum / mass_sum);
                }
                break;
            }
        }

    }

    // return the graph
    return g;
}

// ---------------------------------------------------------
void BCH2D::DrawProfileGraphs()
{
    if (fDrawProfileX) {
        TGraph* graph_profile = CalculateProfileGraph(kProfileX, fProfileXType);
        graph_profile->SetLineColor(fProfileXLineColor);
        graph_profile->SetLineStyle(fProfileXLineStyle);
        graph_profile->Draw("sameL");
        fROOTObjects.push_back(graph_profile);
        std::string xtext = "profile x";
        switch (fProfileXType) {
            case kProfileMode:
                xtext += " (mode)";
                break;
            case kProfileMedian:
                xtext += " (median)";
                break;
            case kProfileMean:
            default:
                xtext += " (mean)";
                break;
        }
        AddLegendEntry(graph_profile, xtext.data(), "L");
    }
    if (fDrawProfileY) {
        TGraph* graph_profile = CalculateProfileGraph(kProfileY, fProfileYType);
        graph_profile->SetLineColor(fProfileYLineColor);
        graph_profile->SetLineStyle(fProfileYLineStyle);
        graph_profile->Draw("sameL");
        fROOTObjects.push_back(graph_profile);
        std::string ytext = "profile y";
        switch (fProfileYType) {
            case kProfileMode:
                ytext += " (mode)";
                break;
            case kProfileMedian:
                ytext += " (median)";
                break;
            case kProfileMean:
            default:
                ytext += " (mean)";
                break;
        }
        AddLegendEntry(graph_profile, ytext.data(), "L");
    }
}

