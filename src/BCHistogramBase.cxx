/*
 * Copyright (C) 2007-2013, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

// ---------------------------------------------------------

#include "BCHistogramBase.h"

#include "BCLog.h"
#include "BCMath.h"

#include <TH1D.h>
#include <TAxis.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TH2.h>
#include <TLine.h>
#include <TPolyLine.h>
#include <TPaveLabel.h>
#include <TLatex.h>
#include <TError.h>
#include <TCanvas.h>
#include <TMarker.h>
#include <TArrow.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TString.h>

#include <math.h>

// ---------------------------------------------------------
BCHistogramBase::BCHistogramBase(TH1 * hist, int dimension)
  : fHistogram(0)
	, fLegend(new TLegend)
	, fBandFillStyle(1001)
	, fMarkerColor(kBlack)
	, fMarkerScale(2)
	, fLogx(false)
	, fLogy(false)
	, fNBands(3)
	, fNSmooth(0)
	, fDrawGlobalMode(true)
	, fDrawGlobalModeArrows(true)
	, fDrawLocalMode(false)
	, fDrawLocalModeArrows(true)
	, fDrawMean(true)
	, fDrawStandardDeviation(true)
	, fDrawLegend(true)
	, fDrawStats(false)
	, fDimension(dimension)
{
	SetHistogram(hist);
	SetColorScheme(kGreenYellowRed);

	fLegend -> SetBorderSize(0);
	fLegend -> SetFillColor(kWhite);
	fLegend -> SetTextAlign(12);
	fLegend -> SetTextFont(62);
	fLegend -> SetTextSize(0.03);
	fROOTObjects.push_back(fLegend);
}

// ---------------------------------------------------------
BCHistogramBase::BCHistogramBase(const BCHistogramBase & other)
  : fHistogram(0)
	, fLegend(new TLegend)
	, fGlobalMode(other.fGlobalMode)
	, fDimension(other.fDimension)
{
	SetHistogram(other.fHistogram);
	
	fLegend -> SetBorderSize(0);
	fLegend -> SetFillColor(kWhite);
	fLegend -> SetTextAlign(12);
	fLegend -> SetTextFont(62);
	fLegend -> SetTextSize(0.03);
	fROOTObjects.push_back(fLegend);

	CopyOptions(other);
}

// ---------------------------------------------------------
void BCHistogramBase::CopyOptions(const BCHistogramBase & other) {
	fBandFillStyle = other.fBandFillStyle;
	fMarkerColor = other.fMarkerColor;
	fMarkerScale = other.fMarkerScale;
	fLogx = other.fLogx;
	fLogy = other.fLogy;
	fNBands = other.fNBands;
	fNSmooth = other.fNSmooth;
	fDrawGlobalMode = other.fDrawGlobalMode;
	fDrawGlobalModeArrows = other.fDrawGlobalModeArrows;
	fDrawLocalMode = other.fDrawLocalMode;
	fDrawLocalModeArrows = other.fDrawLocalModeArrows;
	fDrawMean = other.fDrawMean;
	fDrawStandardDeviation = other.fDrawStandardDeviation;
	fDrawLegend = other.fDrawLegend;
	fDrawStats = other.fDrawStats;
	fBandColors = other.fBandColors;
	fIntervals = other.fIntervals;
}

// ---------------------------------------------------------
BCHistogramBase::~BCHistogramBase()
{
	if (fHistogram)
		delete fHistogram;
   for (unsigned i = 0; i < fROOTObjects.size(); ++i)
      delete fROOTObjects[i];
}

// ---------------------------------------------------------
void BCHistogramBase::SetHistogram(TH1 * hist) {
	if (fHistogram)
		delete fHistogram;

	if (!hist or (fDimension>=0 and hist->GetDimension()!=fDimension)) {
		fHistogram = 0;
		fLocalMode.clear();
		return;
	}

	fHistogram = (TH1*) hist -> Clone(TString::Format("%s_bch",hist->GetName()));
	fDimension = fHistogram->GetDimension();

	// normalize
	double integral = GetHistogram() -> Integral("width");
	if (integral != 0)
		GetHistogram() -> Scale(1./integral);

	// Get local mode
	int b = GetHistogram() -> GetMaximumBin();
	int bx, by, bz;
	GetHistogram() -> GetBinXYZ(b,bx,by,bz);
	fLocalMode.assign(1,GetHistogram()->GetXaxis()->GetBinCenter(bx));
	if (GetHistogram()->GetDimension()>1)
		fLocalMode.push_back(GetHistogram()->GetYaxis()->GetBinCenter(by));
	if (GetHistogram()->GetDimension()>2)
		fLocalMode.push_back(GetHistogram()->GetZaxis()->GetBinCenter(bz));

	// Set Y title, if 1D
	if (GetHistogram()->GetDimension()==1 and strlen(GetHistogram()->GetYaxis()->GetTitle())==0)
		GetHistogram() -> SetYTitle(TString::Format("P(%s|Data)",GetHistogram()->GetXaxis()->GetTitle()));
}

// ---------------------------------------------------------
void BCHistogramBase::SetColorScheme(BCHColorScheme scheme) {
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
		break;
		
	case kBlueOrange:
		AddBandColor(kBlue);
		AddBandColor(kBlue-3);
		AddBandColor(kBlue-1);
		AddBandColor(kBlue-6);
		AddBandColor(kBlue-8);
		AddBandColor(kBlue-9);
		AddBandColor(kBlue-10);
		SetMarkerColor(kOrange);
		break;

	case kRedGreen:
		AddBandColor(kRed);
		AddBandColor(kRed-3);
		AddBandColor(kRed-1);
		AddBandColor(kRed-6);
		AddBandColor(kRed-8);
		AddBandColor(kRed-9);
		AddBandColor(kRed-10);
		SetMarkerColor(kGreen);
		break;

	case kGreenYellowRed:
	default:
		AddBandColor(kGreen);
		AddBandColor(kYellow);
		AddBandColor(kRed);
		AddBandColor(kRed-3);
		AddBandColor(kRed-1);
		AddBandColor(kRed-6);
		SetMarkerColor(kBlack);
		break;

	}

}

// ---------------------------------------------------------
void BCHistogramBase::Smooth(unsigned n) {
	if (n==0)
		return;
	GetHistogram() -> Smooth(fNSmooth);
	double integral = GetHistogram() -> Integral("width");
	if ( integral != 0 )
		GetHistogram() -> Scale(1./integral);
}

// ---------------------------------------------------------
void BCHistogramBase::DrawGlobalMode() {
	double ymin = gPad -> GetUymin();
	double ymax = gPad -> GetUymax();
	double y = ymin + 0.5*(ymax+ymin);
	if (gPad->GetLogy()) {
		ymin = pow(10,ymin);
		ymax = pow(10,ymax);
		y = ymin * pow(ymax/ymin,0.5);
	}
	if (GetHistogram()->GetDimension()>1 and fGlobalMode.size()>1)
		y = fGlobalMode[1];

	if (fDrawGlobalMode and !fGlobalMode.empty()) {
		TMarker * marker_mode = new TMarker(fGlobalMode[0], y, 24);
		marker_mode -> SetMarkerColor(GetMarkerColor());
		marker_mode -> SetMarkerSize(fMarkerScale*gPad->GetWNDC());
		marker_mode -> Draw();
		fROOTObjects.push_back(marker_mode);

		TLegendEntry * le = fLegend -> AddEntry(marker_mode, "global mode", "P");
		le -> SetMarkerStyle(marker_mode->GetMarkerStyle());
		le -> SetMarkerSize(marker_mode->GetMarkerSize());
		le -> SetMarkerColor(marker_mode->GetMarkerColor());

		if (fDrawGlobalModeArrows) {
			TArrow* arrow_mode = new TArrow(marker_mode->GetX(), (gPad->GetLogy() ? marker_mode->GetY()*pow(ymax/ymin,-1.5e-2) : marker_mode->GetY()+(ymax-ymin)*-1.5e-2),
																			marker_mode->GetX(), (gPad->GetLogy() ? ymin*pow(ymax/ymin,3e-2) : ymin+(ymax-ymin)*3e-2),
																			2e-2*gPad->GetWNDC(), "|>");
			arrow_mode -> SetLineColor(marker_mode->GetMarkerColor());
			arrow_mode -> SetFillColor(marker_mode->GetMarkerColor());
			arrow_mode -> Draw();
			fROOTObjects.push_back(arrow_mode);
			
			if (GetHistogram()->GetDimension()>1 and fGlobalMode.size()>1) {
				double xmin = gPad -> GetUxmin();
				double xmax = gPad -> GetUxmax();
				if (gPad->GetLogx()) {
					ymin = pow(10,xmin);
					ymax = pow(10,xmax);
				}
				TArrow* arrow_mode2 = new TArrow((gPad->GetLogx() ? marker_mode->GetX()*pow(xmax/xmin,-1.5e-2) : marker_mode->GetX()+(xmax-xmin)*-1.5e-2), marker_mode->GetY(),
																				 (gPad->GetLogx() ? xmin*pow(xmax/xmin,3e-2) : xmin+(xmax-xmin)*3e-2), marker_mode->GetY(),
																				 2e-2*gPad->GetWNDC(), "|>");
				arrow_mode2 -> SetLineColor(marker_mode->GetMarkerColor());
				arrow_mode2 -> SetFillColor(marker_mode->GetMarkerColor());
				arrow_mode2 -> Draw();
				fROOTObjects.push_back(arrow_mode2);
			}
		}

	}
}

// ---------------------------------------------------------
void BCHistogramBase::DrawLocalMode() {
	double ymin = gPad -> GetUymin();
	double ymax = gPad -> GetUymax();
	double y = ymin + 0.25*(ymax+ymin);
	if (gPad->GetLogy()) {
		ymin = pow(10,ymin);
		ymax = pow(10,ymax);
		y = ymin * pow(ymax/ymin,0.25);
	}
	if (GetHistogram()->GetDimension()>1 and fGlobalMode.size()>1)
		y = fLocalMode[1];

	if (fDrawLocalMode and !fLocalMode.empty()) {
		TMarker * marker_mode = new TMarker(fLocalMode[0], y, 25);
		marker_mode -> SetMarkerColor(GetMarkerColor());
		marker_mode -> SetMarkerSize(fMarkerScale*gPad->GetWNDC());
		marker_mode -> Draw();
		fROOTObjects.push_back(marker_mode);

		TLegendEntry * le = fLegend -> AddEntry(marker_mode, "local mode", "P");
		le -> SetMarkerStyle(marker_mode->GetMarkerStyle());
		le -> SetMarkerSize(marker_mode->GetMarkerSize());
		le -> SetMarkerColor(marker_mode->GetMarkerColor());

		if (fDrawLocalModeArrows) {
			TArrow* arrow_mode = new TArrow(marker_mode->GetX(), (gPad->GetLogy() ? marker_mode->GetY()*pow(ymax/ymin,-1.5e-2) : marker_mode->GetY()+(ymax-ymin)*-1.5e-2),
																			marker_mode->GetX(), (gPad->GetLogy() ? ymin*pow(ymax/ymin,3e-2) : ymin+(ymax-ymin)*3e-2),
																			2e-2*gPad->GetWNDC(), "|>");
			arrow_mode -> SetLineColor(marker_mode->GetMarkerColor());
			arrow_mode -> SetFillColor(marker_mode->GetMarkerColor());
			arrow_mode -> Draw();
			fROOTObjects.push_back(arrow_mode);
			
			if (GetHistogram()->GetDimension()>1 and fGlobalMode.size()>1) {
				double xmin = gPad -> GetUxmin();
				double xmax = gPad -> GetUxmax();
				if (gPad->GetLogx()) {
					ymin = pow(10,xmin);
					ymax = pow(10,xmax);
				}
				TArrow* arrow_mode2 = new TArrow((gPad->GetLogx() ? marker_mode->GetX()*pow(xmax/xmin,-1.5e-2) : marker_mode->GetX()+(xmax-xmin)*-1.5e-2), marker_mode->GetY(),
																				 (gPad->GetLogx() ? xmin*pow(xmax/xmin,3e-2) : xmin+(xmax-xmin)*3e-2), marker_mode->GetY(),
																				 2e-2*gPad->GetWNDC(), "|>");
				arrow_mode2 -> SetLineColor(marker_mode->GetMarkerColor());
				arrow_mode2 -> SetFillColor(marker_mode->GetMarkerColor());
				arrow_mode2 -> Draw();
				fROOTObjects.push_back(arrow_mode2);
			}
		}

	}
	
}

// ---------------------------------------------------------
void BCHistogramBase::DrawMean() {
	double ymin = gPad -> GetUymin();
	double ymax = gPad -> GetUymax();
	double y = ymin + 0.6*(ymax-ymin);
	if (gPad->GetLogy()) {
		ymin = pow(10,ymin);
		ymax = pow(10,ymax);
		y = ymin * pow(ymax/ymin,60e-2);
	}
	if (GetHistogram()->GetDimension()>1)
		y = GetHistogram()->GetMean(2);
	
	if ( fDrawMean ) {
		TMarker * marker_mean = new TMarker(GetHistogram()->GetMean(1), y, 20);
		marker_mean -> SetMarkerColor(GetMarkerColor());
		marker_mean -> SetMarkerSize(fMarkerScale*gPad->GetWNDC());
		marker_mean -> Draw();
		fROOTObjects.push_back(marker_mean);

		TLegendEntry * le = 0;
		if ( fDrawStandardDeviation ) {
			TArrow* arrow_std = new TArrow(marker_mean->GetX()-GetHistogram()->GetRMS(1), marker_mean->GetY(),
																		 marker_mean->GetX()+GetHistogram()->GetRMS(1), marker_mean->GetY(),
																		 0.02*gPad->GetWNDC(), "<|>");
			arrow_std -> SetLineColor(marker_mean->GetMarkerColor());
			arrow_std -> SetFillColor(marker_mean->GetMarkerColor());
			arrow_std -> Draw();
			le = fLegend -> AddEntry(arrow_std, "mean and standard deviation", "PL");
			le -> SetLineColor(arrow_std->GetLineColor());
			
			if (GetHistogram()->GetDimension()>1) {
				TArrow* arrow_std2 = new TArrow(marker_mean->GetX(), marker_mean->GetY()-GetHistogram()->GetRMS(2),
																				marker_mean->GetX(), marker_mean->GetY()+GetHistogram()->GetRMS(2),
																				0.02*gPad->GetWNDC(), "<|>");
				arrow_std2 -> SetLineColor(marker_mean->GetMarkerColor());
				arrow_std2 -> SetFillColor(marker_mean->GetMarkerColor());
				arrow_std2 -> Draw();
			}
		} else
			le = fLegend -> AddEntry(marker_mean, "mean", "P");
		le -> SetMarkerStyle(marker_mean->GetMarkerStyle());
		le -> SetMarkerSize(marker_mean->GetMarkerSize());
		le -> SetMarkerColor(marker_mean->GetMarkerColor());
	}
}

// ---------------------------------------------------------
void BCHistogramBase::DrawLegend() {
	double ymin = gPad -> GetUymin();
	double ymax = gPad -> GetUymax();
	if (gPad->GetLogy()) {
		ymin = pow(10,ymin);
		ymax = pow(10,ymax);
	}
	
	if (fDrawLegend) {
		fHistogram->GetYaxis()->SetRangeUser(ymin, ymax*(1.15+fLegend->GetTextSize()*fLegend->GetNRows())*1.05);
		
		gPad->SetTopMargin(0.02);
		
		fLegend -> SetX1NDC(gPad->GetLeftMargin() + 0.10 * (1 - gPad->GetRightMargin() - gPad->GetLeftMargin()));
		fLegend -> SetX2NDC(1 - gPad->GetRightMargin());
		fLegend -> SetY1NDC(1 - gPad->GetTopMargin() - fLegend->GetTextSize()*fLegend->GetNRows());
		fLegend -> SetY2NDC(1 - gPad->GetTopMargin());
		fLegend -> Draw();
		
		// rescale top margin
		gPad->SetTopMargin(1-fLegend->GetY1NDC()+0.01);
	} else
	 	fHistogram->GetYaxis()->SetRangeUser(ymin, ymax*1.155);
}
