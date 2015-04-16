/*
 * Copyright (C) 2007-2014, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

// ---------------------------------------------------------

#include "BCAux.h"

#include <TStyle.h>
#include <TGaxis.h>
#include <TH2C.h>
#include <TH2S.h>
#include <TH2I.h>
#include <TH2F.h>
#include <TH2D.h>

#include <algorithm>
#include <cmath>
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
TH2* BCAux::Transpose(const TH2* const h, std::string name)
{

    if (h == NULL)
        return NULL;

    if (name.empty())
        name = std::string(h->GetName()) + "_tr";

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
        return new TH2C(name.data(), title.data(), nbins_x, xmin, xmax, nbins_y, ymin, ymax);
    if (dynamic_cast<const TH2S*>(h) != NULL)
        return new TH2S(name.data(), title.data(), nbins_x, xmin, xmax, nbins_y, ymin, ymax);
    if (dynamic_cast<const TH2I*>(h) != NULL)
        return new TH2I(name.data(), title.data(), nbins_x, xmin, xmax, nbins_y, ymin, ymax);
    if (dynamic_cast<const TH2F*>(h) != NULL)
        return new TH2F(name.data(), title.data(), nbins_x, xmin, xmax, nbins_y, ymin, ymax);
    if (dynamic_cast<const TH2D*>(h) != NULL)
        return new TH2D(name.data(), title.data(), nbins_x, xmin, xmax, nbins_y, ymin, ymax);
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
std::string BCAux::SafeName(std::string name)
{
    name.erase(std::remove_if(name.begin(), name.end(), ::isspace), name.end());
    return name;
}

