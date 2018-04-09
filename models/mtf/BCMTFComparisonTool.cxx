/*
 * Copyright (C) 2007-2018, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

// ---------------------------------------------------------

#include "BCMTFComparisonTool.h"

#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TGraphAsymmErrors.h>
#include <TLatex.h>

#include <iostream>

// ---------------------------------------------------------
BCMTFComparisonTool::BCMTFComparisonTool(const std::string& name)
{
    SetName(name);
}

// ---------------------------------------------------------
BCMTFComparisonTool::~BCMTFComparisonTool()
{
    for (int i = 0; i < GetNContributions(); ++i)
        delete fHistogramContainer[i];
}

// ---------------------------------------------------------
void BCMTFComparisonTool::AddContribution(const std::string& name, TH1D hist)
{
    // add name to container
    fNameContainer.push_back(name);

    // add histogram to container
    TH1D* h;
    {
        BCAux::RootSideEffectGuard g;
        h = new TH1D(hist);
    }
    fHistogramContainer.push_back(h);

    // add central value to container
    fCentralValueContainer.push_back(hist.GetMean());

    // add uncertainty to container
    fUncertaintyContainer.push_back(hist.GetRMS());
}

// ---------------------------------------------------------
void BCMTFComparisonTool::AddContribution(const std::string& name, double centralvalue, double uncertainty)
{
    // add name to container
    fNameContainer.push_back(name);

    // add 0 to container
    fHistogramContainer.push_back(NULL);

    // add central value to container
    fCentralValueContainer.push_back(centralvalue);

    // add uncertainty to container
    fUncertaintyContainer.push_back(uncertainty);
}

// ---------------------------------------------------------
void BCMTFComparisonTool::PrintHistograms(const std::string& filename)
{
    // create canvas
    TCanvas c1;
    c1.cd();

    // loop over all histograms
    for (unsigned i = 0; i < fHistogramContainer.size(); ++i) {
        // get histogram
        TH1D* hist = fHistogramContainer.at(i);

        // set color
        hist->SetLineColor(2 + i);

        // draw histogram
        if (i == 0) {
            hist->Draw("HIST");
            std::cout << " here as well." << std::endl;
        } else {
            hist->Draw("SAMEHIST");
            std::cout << " here as well 2." << std::endl;
        }
    }

    // print canvas
    c1.Print(filename.data());
}

// ---------------------------------------------------------
void BCMTFComparisonTool::DrawOverview()
{
    // get number of contributions
    int ncontributions = GetNContributions();

    // create graph
    TGraphAsymmErrors* graph_contributions = new TGraphAsymmErrors(ncontributions);
    fTrash.Put(graph_contributions);
    graph_contributions->SetMarkerStyle(20);
    graph_contributions->SetMarkerSize(1);

    // coordinate system
    double xmin = 0.0;
    double xmax = 0.0;
    double xwidth = 0.0;
    double ymin = -0.5;
    double ymax = double(ncontributions) - 0.5;

    // ---- fill the graph ---- //

    // loop over all contributions
    for (int i = 0; i < ncontributions; ++i) {

        // get summary information
        double centralvalue = fCentralValueContainer.at(i);
        double uncertainty  = fUncertaintyContainer.at(i);

        // update coordinate system
        if ((centralvalue - uncertainty) < xmin || i == 0)
            xmin = centralvalue - uncertainty;
        if ((centralvalue + uncertainty) > xmax || i == 0)
            xmax = centralvalue + uncertainty;
        xwidth = xmax - xmin;

        // set point and error
        graph_contributions->SetPoint(i, centralvalue, double(ncontributions - i - 1));
        graph_contributions->SetPointError(i, uncertainty, uncertainty, 0, 0);
    }

    // ---- do the plotting ---- //

    // create histogram for axes
    TH2D* hist_axes = new TH2D("", Form(";%s;", GetSafeName().c_str()), 1, xmin - 0.25 * xwidth, xmax + 1.75 * xwidth, ncontributions, ymin, ymax);
    fTrash.Put(hist_axes);
    hist_axes->SetStats(kFALSE);
    hist_axes->GetYaxis()->SetNdivisions(0);
    hist_axes->GetYaxis()->SetTitleOffset(1.0);

    // create latex
    TLatex* latex = new TLatex();
    fTrash.Put(latex);
    latex->SetTextSize(0.04);
    if (ncontributions >= 10)
        latex->SetTextSize(0.02);
    latex->SetTextAlign(12);

    // draw
    hist_axes->Draw();
    graph_contributions->Draw("SAMEPZ");

    // loop over all contributions and draw labels
    for (int i = 0; i < ncontributions; ++i) {
        latex->DrawLatex(xmax + 0.25 * xwidth, double(ncontributions - i - 1), fNameContainer.at(i).c_str());
    }
    hist_axes->Draw("SAMEAXIS");
}

// ---------------------------------------------------------
void BCMTFComparisonTool::PrintOverview(const std::string& filename)
{
    // create canvas
    TCanvas c1;
    c1.cd();

    // draw the overview
    DrawOverview();

    // print to file
    c1.Print(filename.data());
}

// ---------------------------------------------------------
