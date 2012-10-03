/*
 * Copyright (C) 2008-2012, Daniel Kollar and Kevin Kroeninger.
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

// ---------------------------------------------------------

#include "BCAux.h"

#include <TStyle.h>

// ---------------------------------------------------------

void BCAux::SetStyle()
{
   // canvases
   gStyle->SetCanvasColor     (kWhite);
   gStyle->SetCanvasBorderMode(0);
   gStyle->SetCanvasDefH      (700);
   gStyle->SetCanvasDefW      (735);

   // pads
   gStyle->SetPadTopMargin   (0.05);
   gStyle->SetPadBottomMargin(0.11);
   gStyle->SetPadLeftMargin  (0.15);
   gStyle->SetPadRightMargin (0.05);
   gStyle->SetPadBorderMode  (0);

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
//   gStyle->SetLabelSize(0.05,   "X");
   gStyle->SetLabelFont(62,     "X");
//   gStyle->SetLabelOffset(0.02, "X");

//   gStyle->SetLabelSize(0.05,   "Y");
   gStyle->SetLabelFont(62,     "Y");
//   gStyle->SetLabelOffset(0.02, "Y");

   // titles
   gStyle->SetTitleFillColor(kWhite);
   gStyle->SetTitleBorderSize(0);
   gStyle->SetTitleFont      (62, "");
   gStyle->SetTitleOffset    (0.0, "");
   gStyle->SetTitleH         (0.07);

   gStyle->SetTitleFont      (62,   "X");
//   gStyle->SetTitleSize      (0.06, "X");
   gStyle->SetTitleOffset    (1.1,  "X");

   gStyle->SetTitleFont      (62,   "Y");
//   gStyle->SetTitleSize      (0.06, "Y");
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

// ---------------------------------------------------------

