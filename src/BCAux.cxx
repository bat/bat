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
   gStyle->SetCanvasDefW      ((int)(700.*(1.-gStyle->GetPadTopMargin()-gStyle->GetPadBottomMargin())/(1.-gStyle->GetPadLeftMargin()-gStyle->GetPadRightMargin())));

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
   gStyle->SetPalette(1,0);

}
