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
#include "config.h"

#include <TCanvas.h>
#include <TDirectory.h>
#include <TH1.h>
#include <TH2.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TList.h>
#include <TPad.h>

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>

// ---------------------------------------------------------
void BCAux::SetStyle()
{
    BCLog::OutWarning("BCAux::SetStyle() is deprecated and no longer does anything. Please do not use it.");
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
        else if (ext != ".pdf" && ext != ".ps")
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

    TH2* ht = OwnClone(h);
    ht->SetBins(nbins_x, xmin, xmax, nbins_y, ymin, ymax);
    ht->SetNameTitle(newName.data(), title.data());

    return ht;
}

// ---------------------------------------------------------
BCAux::BCRange BCAux::RangeType(double xmin, double xmax)
{
    if (xmin > xmax)
        // return -RangeType(xmax,xmin);
        return BCAux::kReverseRange;
    if (xmin == xmax)
        return BCAux::kEmptyRange;
    if (std::isfinite(xmin) && std::isfinite(xmax))
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
namespace { static const char* ROOToptions = "HISTSAME"; }

// ---------------------------------------------------------
void BCAux::SetKnowledgeUpdateDrawingStyle(BCH1D& prior, BCH1D& posterior, BCAux::BCKnowledgeUpdateDrawingStyle style)
{
    switch (style) {
        case BCAux::kKnowledgeUpdateDetailedPosterior:
            prior.SetDrawGlobalMode(false);
            prior.SetDrawLocalMode(false);
            prior.SetDrawMean(false);
            prior.SetDrawMedian(false);
            prior.SetDrawLegend(false);
            prior.SetBandType(BCH1D::kNoBands);
            prior.SetROOToptions(ROOToptions);
            prior.SetLineColor(13);
            prior.SetMarkerColor(13);
            prior.SetNLegendColumns(1);
            prior.SetColorScheme(BCHistogramBase::kBlackWhite);
            prior.SetLineStyle(1);
            prior.SetLineColor(kGray + 1);
            prior.SetLineWidth(2);
            posterior.CopyOptions(prior);
            posterior.SetDrawGlobalMode(true);
            posterior.SetBandType(BCH1D::kSmallestInterval);
            posterior.SetNBands(3);
            posterior.SetColorScheme(BCHistogramBase::kGreenYellowRed);
            posterior.SetLineWidth(1);
            break;

        case BCAux::kKnowledgeUpdateDetailedPrior:
            SetKnowledgeUpdateDrawingStyle(posterior, prior, kKnowledgeUpdateDetailedPosterior);
            break;

        case BCAux::kKnowledgeUpdateDefaultStyle:
        default:
            posterior.SetDrawGlobalMode(false);
            posterior.SetDrawLocalMode(false);
            posterior.SetDrawMean(false);
            posterior.SetDrawMedian(false);
            posterior.SetDrawLegend(false);
            posterior.SetBandType(BCH1D::kNoBands);
            posterior.SetROOToptions(ROOToptions);
            posterior.SetColorScheme(BCHistogramBase::kBlackWhite);
            posterior.SetNLegendColumns(1);
            posterior.SetLineWidth(2);
            prior.CopyOptions(posterior);
            prior.SetLineColor(kGray + 1);
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
            prior.SetBandFillStyle(-3);
            prior.SetNBands(1);
            prior.SetNSmooth(0);
            prior.SetROOToptions(ROOToptions);
            prior.SetLineColor(13);
            prior.SetMarkerColor(13);
            prior.SetLocalModeMarkerStyle(25);
            prior.SetNLegendColumns(1);
            prior.SetColorScheme(BCHistogramBase::kBlackWhite);
            prior.SetLineStyle(1);
            prior.SetLineWidth(1);
            posterior.CopyOptions(prior);
            posterior.SetNBands(3);
            posterior.SetBandFillStyle(1001);
            posterior.SetColorScheme(BCHistogramBase::kGreenYellowRed);
            posterior.SetLocalModeMarkerStyle(21);
            break;

        case BCAux::kKnowledgeUpdateDetailedPrior:
            SetKnowledgeUpdateDrawingStyle(posterior, prior, kKnowledgeUpdateDetailedPosterior);
            break;

        case BCAux::kKnowledgeUpdateDefaultStyle:
        default:
            posterior.SetDrawGlobalMode(false);
            posterior.SetDrawLocalMode(true, false);
            posterior.SetDrawMean(false);
            posterior.SetDrawLegend(false);
            posterior.SetBandType(BCH2D::kSmallestInterval);
            posterior.SetBandFillStyle(-3);
            posterior.SetNBands(1);
            posterior.SetNSmooth(0);
            posterior.SetROOToptions(ROOToptions);
            posterior.SetColorScheme(BCHistogramBase::kBlackWhite);
            posterior.SetNLegendColumns(1);
            posterior.SetLineWidth(2);
            posterior.SetLocalModeMarkerStyle(21);
            prior.CopyOptions(posterior);
            prior.SetLineColor(kGray + 1);
            prior.SetMarkerColor(kGray + 1);
            break;
    }
}


// ---------------------------------------------------------
void BCAux::DrawKnowledgeUpdate(BCHistogramBase& prior, BCHistogramBase& posterior, bool draw_prior_first, BCTrash<TObject>& trash)
{
    if (prior.GetDimension() != posterior.GetDimension()) {
        BCLog::OutError("BCAux::DrawKnowledgeUpdate : prior and posterior dimension do not match.");
        return;
    }

    if (!prior.Valid() || !posterior.Valid()) {
        if (!prior.Valid() && posterior.Valid())
            BCLog::OutError("BCAux::DrawKnowledgeUpdate : prior invalid.");
        else if (prior.Valid() && !posterior.Valid())
            BCLog::OutError("BCAux::DrawKnowledgeUpdate : posterior invalid.");
        else
            BCLog::OutError("BCAux::DrawKnowledgeUpdate : prior and posterior invalid.");
        return;
    }

    gPad->SetLogx(prior.GetLogx() || posterior.GetLogx());
    gPad->SetLogy(prior.GetLogy() || posterior.GetLogy());
    gPad->SetGridx(prior.GetGridx() || posterior.GetGridx());
    gPad->SetGridx(prior.GetGridy() || posterior.GetGridy());
    if (prior.GetDimension() > 2)
        gPad->SetLogy(prior.GetLogz() || posterior.GetLogz());

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
    TH2D* h2_axes;
    {
        RootSideEffectGuard g;
        h2_axes = new TH2D(Form("h2_axes_knowledge_update_%s_%s", prior.GetHistogram()->GetName(), posterior.GetHistogram()->GetName()),
                           Form(";%s;%s", prior.GetHistogram()->GetXaxis()->GetTitle(), prior.GetHistogram()->GetYaxis()->GetTitle()),
                           10, minx, maxx, 10, miny, maxy);
    }
    trash.Put(h2_axes);
    h2_axes->SetStats(false);
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

    if (prior.GetHistogram()->GetDimension() > 1 || prior.GetLegend().GetNRows() > 0)
        prior.GetLegend().SetHeader(prior.GetHistogram()->GetTitle());
    else
        prior.GetLegend().AddEntry(prior.GetHistogram(), 0, "L");
    if (posterior.GetLegend().GetNRows() > 0)
        posterior.GetLegend().SetHeader(posterior.GetHistogram()->GetTitle());
    else
        posterior.GetLegend().AddEntry(posterior.GetHistogram(), 0, "L");


    double y1ndc_prior = prior.ResizeLegend();
    double y1ndc_posterior = posterior.ResizeLegend();

    // Draw prior legend on top left
    // Draw posterior legend on top right
#if ROOTVERSION >= 6000000
    prior.GetLegend().SetX2(prior.GetLegend().GetX1() + 45e-2 * (prior.GetLegend().GetX2() - prior.GetLegend().GetX1()));
    posterior.GetLegend().SetX1(posterior.GetLegend().GetX1() + 55e-2 * (posterior.GetLegend().GetX2() - posterior.GetLegend().GetX1()));
#else
    prior.GetLegend().SetX2NDC(prior.GetLegend().GetX1NDC() + 45e-2 * (prior.GetLegend().GetX2NDC() - prior.GetLegend().GetX1NDC()));
    posterior.GetLegend().SetX1NDC(posterior.GetLegend().GetX1NDC() + 55e-2 * (posterior.GetLegend().GetX2NDC() - posterior.GetLegend().GetX1NDC()));
#endif

    prior.GetLegend().Draw();
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

// ---------------------------------------------------------
BCAux::RootSideEffectGuard::RootSideEffectGuard() :
    fDirectory(gDirectory)
{
    gDirectory = NULL;
}

// ---------------------------------------------------------
BCAux::RootSideEffectGuard::~RootSideEffectGuard()
{
    gDirectory = fDirectory;
}
