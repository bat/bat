/*
 * Copyright (C) 2007-2015, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

// ---------------------------------------------------------

#include "BCAux.h"
#include "BCLog.h"

#include <TCanvas.h>
#include <TGaxis.h>
#include <TH1.h>
#include <TH2C.h>
#include <TH2S.h>
#include <TH2I.h>
#include <TH2F.h>
#include <TH2D.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TList.h>
#include <TPad.h>
#include <TStyle.h>

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>

// ---------------------------------------------------------

void BCAux::SetStyle()
{
    // pads
    gStyle->SetPadTopMargin   (0.05);
    gStyle->SetPadBottomMargin(0.11);
    gStyle->SetPadLeftMargin  (0.15);
    gStyle->SetPadRightMargin (0.10);
    gStyle->SetPadBorderMode  (0);

    // canvases
    gStyle->SetCanvasColor     (kWhite);
    gStyle->SetCanvasBorderMode(0);
    gStyle->SetCanvasDefH      (700);
    gStyle->SetCanvasDefW      ((int)(700.*(1. - gStyle->GetPadTopMargin() - gStyle->GetPadBottomMargin()) / (1. - gStyle->GetPadLeftMargin() - gStyle->GetPadRightMargin())));

    // Frames
    gStyle->SetFrameFillStyle (0);
    gStyle->SetFrameFillColor (kWhite);
    gStyle->SetFrameLineColor (kBlack);
    gStyle->SetFrameLineStyle (0);
    gStyle->SetFrameLineWidth (1);
    gStyle->SetFrameBorderMode(0);

    // histograms
    gStyle->SetHistFillColor(kWhite);
    gStyle->SetHistFillStyle(0);
    gStyle->SetHistLineColor(kBlack);
    gStyle->SetHistLineStyle(0);
    gStyle->SetHistLineWidth(1);
    gStyle->SetStripDecimals(kFALSE);

    // set decimals
    TGaxis::SetMaxDigits(4);

    // options
    gStyle->SetOptTitle(0);

    // lines
    gStyle->SetLineColor(kBlack);
    gStyle->SetLineStyle(1);
    gStyle->SetLineWidth(1);

    // markers
    gStyle->SetMarkerStyle(kFullCircle);
    gStyle->SetMarkerSize (1.0);

    // functions
    gStyle->SetFuncColor(kBlack);
    gStyle->SetFuncStyle(0);
    gStyle->SetFuncWidth(2);

    // labels
    gStyle->SetLabelFont(62,     "X");
    gStyle->SetLabelFont(62,     "Y");

    // titles
    gStyle->SetTitleFillColor(kWhite);
    gStyle->SetTitleBorderSize(0);
    gStyle->SetTitleFont      (62, "");
    gStyle->SetTitleOffset    (0.0, "");
    gStyle->SetTitleH         (0.07);

    gStyle->SetTitleFont      (62,   "X");
    gStyle->SetTitleOffset    (1.1,  "X");

    gStyle->SetTitleFont      (62,   "Y");
    gStyle->SetTitleOffset    (1.8,  "Y");

    // ticks
    gStyle->SetTickLength(0.03);

    // statistics box
    gStyle->SetStatFont (62);
    gStyle->SetStatColor(kWhite);
    gStyle->SetStatH    (0.20);
    gStyle->SetStatW    (0.20);
    gStyle->SetStatX    (0.965);
    gStyle->SetStatY    (0.90);

    // palette
    gStyle->SetPalette(1, 0);

}

// ---------------------------------------------------------
void BCAux::DefaultToPDF(std::string& filename)
{
    if (filename.empty())
        return;

    size_t ext_pos = filename.find_last_of(".");
    if (ext_pos == std::string::npos)
        filename += ".pdf";
    else {
        std::string ext = filename.substr(ext_pos);
        std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);
        if (ext == ".")
            filename += "pdf";
        else if (ext != ".pdf" and ext != ".ps")
            filename += ".pdf";
    }
}

// ---------------------------------------------------------
TH2* BCAux::Transpose(const TH2* const h, const std::string& name)
{

    if (h == NULL)
        return NULL;

    std::string newName(name);
    if (newName.empty())
        newName = std::string(h->GetName()) + "_tr";

    int nbins_x = h->GetNbinsY();
    double xmin = h->GetYaxis()->GetXmin();
    double xmax = h->GetYaxis()->GetXmax();
    std::string xtitle = h->GetYaxis()->GetTitle();

    int nbins_y = h->GetNbinsX();
    double ymin = h->GetXaxis()->GetXmin();
    double ymax = h->GetXaxis()->GetXmax();
    std::string ytitle = h->GetXaxis()->GetTitle();

    std::string title = std::string(h->GetTitle()) + ";" + xtitle + ";" + ytitle + ";" + h->GetZaxis()->GetTitle();

    if (dynamic_cast<const TH2C*>(h) != NULL)
        return new TH2C(newName.data(), title.data(), nbins_x, xmin, xmax, nbins_y, ymin, ymax);
    if (dynamic_cast<const TH2S*>(h) != NULL)
        return new TH2S(newName.data(), title.data(), nbins_x, xmin, xmax, nbins_y, ymin, ymax);
    if (dynamic_cast<const TH2I*>(h) != NULL)
        return new TH2I(newName.data(), title.data(), nbins_x, xmin, xmax, nbins_y, ymin, ymax);
    if (dynamic_cast<const TH2F*>(h) != NULL)
        return new TH2F(newName.data(), title.data(), nbins_x, xmin, xmax, nbins_y, ymin, ymax);
    if (dynamic_cast<const TH2D*>(h) != NULL)
        return new TH2D(newName.data(), title.data(), nbins_x, xmin, xmax, nbins_y, ymin, ymax);
    return NULL;
}

// ---------------------------------------------------------
BCAux::BCRange BCAux::RangeType(double xmin, double xmax)
{
    if (xmin > xmax)
        // return -RangeType(xmax,xmin);
        return BCAux::kReverseRange;
    if (xmin == xmax)
        return BCAux::kEmptyRange;
    if (std::isfinite(xmin) and std::isfinite(xmax))
        return BCAux::kFiniteRange;
    if (std::isfinite(xmax))
        return BCAux::kNegativeInfiniteRange;
    if (std::isfinite(xmin))
        return BCAux::kPositiveInfiniteRange;
    return BCAux::kInfiniteRange;
}

// ---------------------------------------------------------
void BCAux::MakeFinite(double& xmin, double& xmax)
{
    if (!std::isfinite(xmin))
        xmin = -std::numeric_limits<double>::max();
    if (!std::isfinite(xmax))
        xmax = +std::numeric_limits<double>::max();
}

// ---------------------------------------------------------
std::string BCAux::SafeName(const std::string& name)
{
    std::string res(name);
    res.erase(std::remove_if(res.begin(), res.end(), std::not1(std::ptr_fun(BCAux::AllowedCharacter))), res.end());
    return res;
}

// ---------------------------------------------------------
bool BCAux::AllowedCharacter(char c)
{
    if (::isalnum(c))
        return true;
    if (c == '_')
        return true;
    return false;
}

// ---------------------------------------------------------
void BCAux::SetKnowledgeUpdateDrawingStyle(BCH1D& prior, BCH1D& posterior, BCAux::BCKnowledgeUpdateDrawingStyle style)
{
    switch (style) {

        case BCAux::kKnowledgeUpdateDetailedPosterior:
            // 1D
            prior.SetDrawGlobalMode(false);
            prior.SetDrawLocalMode(false);
            prior.SetDrawMean(false);
            prior.SetDrawMedian(false);
            prior.SetDrawLegend(false);
            prior.SetNBands(0);
            prior.SetBandType(BCH1D::kNoBands);
            prior.SetROOToptions("same");
            prior.SetLineColor(13);
            prior.SetMarkerColor(13);
            prior.SetNLegendColumns(1);
            posterior.CopyOptions(prior);
            posterior.SetDrawGlobalMode(true);
            posterior.SetBandType(BCH1D::kSmallestInterval);
            posterior.SetNBands(3);
            posterior.SetNLegendColumns(1);
            posterior.SetColorScheme(BCHistogramBase::kGreenYellowRed);
            posterior.SetLineColor(kBlack);
            posterior.SetMarkerColor(kBlack);
            break;

        case BCAux::kKnowledgeUpdateDetailedPrior:
            posterior.SetDrawGlobalMode(false);
            posterior.SetDrawLocalMode(false);
            posterior.SetDrawMean(false);
            posterior.SetDrawMedian(false);
            posterior.SetDrawLegend(false);
            posterior.SetNBands(0);
            posterior.SetBandType(BCH1D::kNoBands);
            posterior.SetROOToptions("same");
            posterior.SetLineColor(13);
            posterior.SetMarkerColor(13);
            posterior.SetNLegendColumns(1);
            prior.CopyOptions(posterior);
            prior.SetDrawGlobalMode(true);
            prior.SetBandType(BCH1D::kSmallestInterval);
            prior.SetNBands(3);
            prior.SetNLegendColumns(1);
            prior.SetColorScheme(BCHistogramBase::kGreenYellowRed);
            prior.SetLineColor(kBlack);
            prior.SetMarkerColor(kBlack);
            break;

        case BCAux::kKnowledgeUpdateDefaultStyle:
        default:
            // 1D
            prior.SetDrawGlobalMode(false);
            prior.SetDrawLocalMode(false);
            prior.SetDrawMean(false);
            prior.SetDrawMedian(false);
            prior.SetDrawLegend(false);
            prior.SetNBands(0);
            prior.SetBandType(BCH1D::kNoBands);
            prior.SetROOToptions("same");
            prior.SetLineColor(kRed);
            prior.SetMarkerColor(kRed);
            prior.SetNLegendColumns(1);
            posterior.CopyOptions(prior);
            posterior.SetNLegendColumns(1);
            posterior.SetLineColor(kBlue);
            posterior.SetMarkerColor(kBlue);
            break;
    }
}

// ---------------------------------------------------------
void BCAux::SetKnowledgeUpdateDrawingStyle(BCH2D& prior, BCH2D& posterior, BCAux::BCKnowledgeUpdateDrawingStyle style)
{
    switch (style) {

        case BCAux::kKnowledgeUpdateDetailedPosterior:
            prior.SetDrawGlobalMode(false);
            prior.SetDrawLocalMode(true, false);
            prior.SetDrawMean(false);
            prior.SetDrawLegend(false);
            prior.SetBandType(BCH2D::kSmallestInterval);
            prior.SetBandFillStyle(-1);
            prior.SetNBands(1);
            prior.SetNSmooth(0);
            prior.SetROOToptions("same");
            prior.SetLineColor(13);
            prior.SetMarkerColor(13);
            prior.SetLocalModeMarkerStyle(25);
            prior.SetNLegendColumns(1);
            posterior.CopyOptions(prior);
            posterior.SetNBands(3);
            posterior.SetBandFillStyle(1001);
            posterior.SetColorScheme(BCHistogramBase::kGreenYellowRed);
            posterior.SetLineColor(kBlack);
            posterior.SetMarkerColor(kBlack);
            posterior.SetLocalModeMarkerStyle(21);
            break;

        case BCAux::kKnowledgeUpdateDetailedPrior:
            // 2D
            posterior.SetDrawGlobalMode(false);
            posterior.SetDrawLocalMode(true, false);
            posterior.SetDrawMean(false);
            posterior.SetDrawLegend(false);
            posterior.SetBandType(BCH2D::kSmallestInterval);
            posterior.SetBandFillStyle(-1);
            posterior.SetNBands(1);
            posterior.SetNSmooth(0);
            posterior.SetROOToptions("same");
            posterior.SetLineColor(13);
            posterior.SetMarkerColor(13);
            posterior.SetLocalModeMarkerStyle(25);
            posterior.SetNLegendColumns(1);
            prior.CopyOptions(posterior);
            prior.SetNBands(3);
            prior.SetColorScheme(BCHistogramBase::kGreenYellowRed);
            prior.SetBandFillStyle(1001);
            prior.SetLineColor(kBlack);
            prior.SetMarkerColor(kBlack);
            prior.SetLocalModeMarkerStyle(21);
            break;

        case BCAux::kKnowledgeUpdateDefaultStyle:
        default:
            // 2D
            prior.SetDrawGlobalMode(false);
            prior.SetDrawLocalMode(true, false);
            prior.SetDrawMean(false);
            prior.SetDrawLegend(false);
            prior.SetBandType(BCH2D::kSmallestInterval);
            prior.SetBandFillStyle(-1);
            prior.SetNBands(1);
            prior.SetNSmooth(0);
            prior.SetROOToptions("same");
            prior.SetLineColor(kRed);
            prior.SetMarkerColor(kRed);
            prior.SetNLegendColumns(1);
            posterior.CopyOptions(prior);
            posterior.SetLineColor(kBlue);
            posterior.SetMarkerColor(kBlue);
            break;
    }
}


// ---------------------------------------------------------
void BCAux::DrawKnowledgeUpdate(BCHistogramBase& prior, BCHistogramBase& posterior, bool draw_prior_first)
{
    if (prior.GetDimension() != posterior.GetDimension()) {
        BCLog::OutError("BCAux::DrawKnowledgeUpdate : prior and posterior dimension do not match.");
        return;
    }

    if (!prior.Valid() or !posterior.Valid()) {
        if (!prior.Valid() and posterior.Valid())
            BCLog::OutError("BCAux::DrawKnowledgeUpdate : prior invalid.");
        else if (prior.Valid() and !posterior.Valid())
            BCLog::OutError("BCAux::DrawKnowledgeUpdate : posterior invalid.");
        else
            BCLog::OutError("BCAux::DrawKnowledgeUpdate : prior and posterior invalid.");
        return;
    }

    gPad->SetLogx(prior.GetLogx() or posterior.GetLogx());
    gPad->SetLogy(prior.GetLogy() or posterior.GetLogy());
    gPad->SetGridx(prior.GetGridx() or posterior.GetGridx());
    gPad->SetGridx(prior.GetGridy() or posterior.GetGridy());
    if (prior.GetDimension() > 2)
        gPad->SetLogy(prior.GetLogz() or posterior.GetLogz());

    // get ranges
    double minx = std::min<double>(prior.GetHistogram()->GetXaxis()->GetXmin(), posterior.GetHistogram()->GetXaxis()->GetXmin());
    double maxx = std::max<double>(prior.GetHistogram()->GetXaxis()->GetXmax(), posterior.GetHistogram()->GetXaxis()->GetXmax());

    double miny = 0;
    double maxy = 0;

    if (prior.GetDimension() == 1) {
        miny = 0.0;
        maxy = 1.1 * std::max<double>(prior.GetHistogram()->GetMaximum(), posterior.GetHistogram()->GetMaximum());
        if (gPad->GetLogy()) {
            miny = 0.5 * std::min<double>(prior.GetHistogram()->GetMinimum(0), posterior.GetHistogram()->GetMinimum(0));
            maxy *= 2;
        }
    } else {
        miny = std::min<double>(prior.GetHistogram()->GetYaxis()->GetXmin(), posterior.GetHistogram()->GetYaxis()->GetXmin());
        maxy = std::max<double>(prior.GetHistogram()->GetYaxis()->GetXmax(), posterior.GetHistogram()->GetYaxis()->GetXmax());
    }

    // draw axes
    TH2D* h2_axes = new TH2D(Form("h2_axes_knowledge_update_%s_%s", prior.GetHistogram()->GetName(), posterior.GetHistogram()->GetName()),
                             Form(";%s;%s", prior.GetHistogram()->GetXaxis()->GetTitle(), prior.GetHistogram()->GetYaxis()->GetTitle()),
                             10, minx, maxx, 10, miny, maxy);
    h2_axes->SetStats(false);
    h2_axes->GetXaxis()->SetNdivisions(508);
    if (prior.GetDimension() > 1)
        h2_axes->GetYaxis()->SetNdivisions(508);
    h2_axes->Draw();

    // turn off legends
    prior.SetDrawLegend(false);
    posterior.SetDrawLegend(false);

    // Draw histograms
    // ROOT options for both prior and posterior should contain "same"
    // (as they do by default)
    if (!draw_prior_first)
        posterior.Draw();
    prior.Draw();
    if (draw_prior_first)
        posterior.Draw();

    // create / draw legend(s)

    if (prior.GetLegend().GetNRows() > 0)
        prior.GetLegend().SetHeader(prior.GetHistogram()->GetTitle());
    else
        prior.GetLegend().AddEntry(prior.GetHistogram(), 0, "L");
    if (posterior.GetLegend().GetNRows() > 0)
        posterior.GetLegend().SetHeader(posterior.GetHistogram()->GetTitle());
    else
        posterior.GetLegend().AddEntry(posterior.GetHistogram(), 0, "L");


    // Draw prior legend on top left
    double y1ndc_prior = prior.ResizeLegend();
    prior.GetLegend().SetX2NDC(prior.GetLegend().GetX1NDC() + 45e-2 * (prior.GetLegend().GetX2NDC() - prior.GetLegend().GetX1NDC()));
    prior.GetLegend().Draw();

    // Draw posterior legend on top right
    double y1ndc_posterior = posterior.ResizeLegend();
    posterior.GetLegend().SetX1NDC(posterior.GetLegend().GetX1NDC() + 55e-2 * (posterior.GetLegend().GetX2NDC() - posterior.GetLegend().GetX1NDC()));
    posterior.GetLegend().Draw();

    gPad->SetTopMargin(1 - std::min<double>(y1ndc_prior, y1ndc_posterior) + 0.01);

    gPad->RedrawAxis();
    gPad->Update();
}

// ---------------------------------------------------------
unsigned BCAux::PrintPlots(std::vector<BCH1D>& h1, std::vector<BCH2D>& h2, const std::string& filename, unsigned hdiv, unsigned vdiv)
{
    const unsigned nplots = h1.size() + h2.size();
    if (nplots == 0) {
        BCLog::OutWarning("BCAux::PrintPlots : No plots to print");
        return 0;
    }

    BCLog::OutSummary(Form("Printing all marginalized distributions (%lu x 1D + %lu x 2D = %u) into file %s", h1.size(), h2.size(), nplots, filename.c_str()));
    // give out warning if too many plots
    if (nplots > 100)
        BCLog::OutDetail("This can take a while...");

    // setup the canvas and file
    if (hdiv < 1) hdiv = 1;
    if (vdiv < 1) vdiv = 1;

    int c_width  = 297 * 4;
    int c_height = 210 * 4;
    if (hdiv < vdiv)
        std::swap(c_width, c_height);

    TCanvas c("c", "canvas", c_width, c_height);
    c.Divide(hdiv, vdiv);

    // open file
    c.Print((filename + "[").data());

    unsigned n = 0;

    // 1D
    for (unsigned i = 0; i < h1.size(); ++i) {
        // if current page is full, switch to new page
        if (i != 0 && i % (hdiv * vdiv) == 0) {
            c.Print(filename.c_str());
            c.Clear("D");
        }

        // go to next pad
        c.cd(i % (hdiv * vdiv) + 1)->ResetAttPad();

        h1[i].Draw();

        if (++n % 100 == 0)
            BCLog::OutDetail(Form(" --> %d plots done", n));
    }
    if (h1.size() > 0) {
        c.Print(filename.c_str());
        c.Clear("D");
    }

    // 2D
    for (unsigned i = 0; i < h2.size(); ++i) {
        // if current page is full, switch to new page
        if (i != 0 && i % (hdiv * vdiv) == 0) {
            c.Print(filename.c_str());
            c.Clear("D");
        }

        // go to next pad
        c.cd(i % (hdiv * vdiv) + 1)->ResetAttPad();

        h2[i].Draw();

        if (++n % 100 == 0)
            BCLog::OutDetail(Form(" --> %d plots done", n));
    }
    if (h2.size() > 0) {
        c.Print(filename.c_str());
        c.Clear("D");
    }

    // close file
    c.Print(std::string( filename + "]").c_str());

    if (nplots > 100 && nplots % 100 != 0)
        BCLog::OutDetail(Form(" --> %d plots done", nplots));

    // return total number of drawn histograms
    return nplots;
}
