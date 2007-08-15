#include "BCH1D.h"

#include <TCanvas.h>
#include <TLine.h>
#include <TMarker.h>
#include <TLegend.h>

// ---------------------------------------------------------

BCH1D::BCH1D()
{
	fHistogram = 0;
	fDefaultCLLimit = 95.; // %

	// suppress the ROOT Info printouts
	gErrorIgnoreLevel=1;
}

// --------------------------------------------------------- 

BCH1D::~BCH1D()
{
	if (fHistogram)
		delete fHistogram;
}

// --------------------------------------------------------- 

double BCH1D::GetMode()
{
	return fHistogram -> GetBinCenter(fHistogram -> GetMaximumBin()); 
}

// --------------------------------------------------------- 

double BCH1D::GetQuantile(double probabilitysum)
{

  int nquantiles = 1; 
  double quantiles[1]; 
  double probsum[1]; 

  probsum[0] = probabilitysum; 

  fHistogram -> GetQuantiles(nquantiles, quantiles, probsum); 

  return quantiles[0]; 

}

// --------------------------------------------------------- 

double BCH1D::GetIntegral(double valuemin, double valuemax) 
{

  double integral = 0; 

  int binmin = fHistogram -> FindBin(valuemin); 
  int binmax = fHistogram -> FindBin(valuemax); 

  integral = fHistogram -> Integral(binmin, binmax); 

  return integral; 

}

// --------------------------------------------------------- 

double BCH1D::GetPValue(double probability)
{
  BCLog::Out(BCLog::warning,BCLog::warning, "BCH1D::GetPValue.  Kevin !!! shouldn't we change this valuemax from bin center to upper edge?");
  double valuemax = fHistogram -> GetBinCenter(fHistogram -> GetNbinsX()); 
  
  double integral = this -> GetIntegral(probability, valuemax); 

  return integral; 

}

// --------------------------------------------------------- 

void BCH1D::SetDefaultCLLimit(double limit)
{
	if(limit>=100. || limit<68.)
		BCLog::Out(BCLog::warning,BCLog::warning,
			Form("BCH1D::SetDefaultCLLimit. Value %.1f out of range. Keeping default %.1f%% CL limit.",limit,fDefaultCLLimit));
	else
		fDefaultCLLimit = limit;
}

// --------------------------------------------------------- 

void BCH1D::Print(char * filename, int options, double ovalue)
{
	// create temporary canvas
	TCanvas * canvas = new TCanvas();
	canvas -> cd();

	double min, max;
	double mode = this->GetMode();

	TLine * line;

	switch(options)
	{
		// standard 68% band drawing
		// if mode is outside the band, draw limit
		case 0:
			if (TMath::Abs(ovalue) >= 100 || TMath::Abs(ovalue) < 68)
			{
				min = this -> GetQuantile(.16);
				max = this -> GetQuantile(.84);
				if (mode > max)
				{
					min = this -> GetQuantile(1.-(double)fDefaultCLLimit/100.);
					max = fHistogram->GetXaxis()->GetXmax();
					ovalue = fDefaultCLLimit;
				}
				else if (mode < min)
				{
					min = fHistogram->GetXaxis()->GetXmin();
					max = this -> GetQuantile((double)fDefaultCLLimit/100.);
					ovalue = -fDefaultCLLimit;
				}
			}
			else if(ovalue < 0)
			{
				min = fHistogram->GetXaxis()->GetXmin();
				max = this -> GetQuantile(-ovalue/100.);
			}
			else
			{
				min = this -> GetQuantile((1-ovalue)/100.);
				max = fHistogram->GetXaxis()->GetXmax();
			}

			// do the drawing
			this->PrintShadedLimits(mode, min, max, ovalue);

	      break;

		// line at value
		case 1:
			fHistogram ->Draw();
			line = new TLine();
			line -> SetLineColor(kRed);
			line -> DrawLine(ovalue, 0., ovalue, 1.05 * fHistogram -> GetMaximum());

	      break;

		// smallest interval
		case 2:
			if(ovalue<50) // just to ensure there's some sense in the number
				ovalue = 68.; // default is 68%
			this->GetSmallestInterval(min, max, ovalue/100.);
			this->PrintShadedLimits(mode, min, max, 0);

	      break;

		// bad options
		default:
			BCLog::Out(BCLog::warning, BCLog::warning, Form("BCH1D::Print. Invalid option %d",options));
			break;
	}

	// print to file 
	canvas -> Print(filename); 
}

// --------------------------------------------------------- 

void BCH1D::PrintShadedLimits(double mode, double min, double max, double limit)
{
	double maximum = fHistogram -> GetMaximum();

	double x0 = mode;
	double y0 = fHistogram->GetBinContent( fHistogram->FindBin(mode) );

	double ysize = maximum*1.2;

	double xmin=fHistogram->GetXaxis()->GetXmin();
	double xmax=fHistogram->GetXaxis()->GetXmax();

	// draw histogram with axes first
	TH2D * hsc = new TH2D(Form("h2_%s",fHistogram->GetName()),"",
			10,xmin,xmax,10,0.,ysize);
	hsc -> SetStats(0);
	hsc -> SetXTitle(fHistogram->GetXaxis()->GetTitle());
	hsc -> SetYTitle(fHistogram->GetYaxis()->GetTitle());
	hsc -> Draw();

	// draw histogram
	fHistogram -> Draw("same");

	// draw yellow shaded region between min and max
	TH1D * hist_shaded =  (TH1D *) fHistogram -> Clone();
	hist_shaded -> SetFillStyle(1001);
	hist_shaded -> SetFillColor(kYellow);

	// remove entries on the sides
	int nbins = hist_shaded -> GetNbinsX();

	int lowerborder = hist_shaded -> FindBin(min);
	for (int i = 1; i < lowerborder; i++)
		hist_shaded -> SetBinContent(i, 0.);

	int upperborder = hist_shaded -> FindBin(max) + 1;
	for (int i = upperborder; i <= nbins; i++)
		hist_shaded -> SetBinContent(i, 0.);

	// draw shaded histogram
	hist_shaded -> Draw("same");

	gPad->RedrawAxis();

	// draw line and triangle for mode
	TLine * line;
	TPolyLine * tmax;
	if(TMath::Abs(limit)<50) // just to ensure there's some sense in the number
	{
		line = new TLine();
		line -> SetLineStyle(2);
		line -> DrawLine(x0, 0., x0, y0);
	}

	double dx = 0.01*(xmax-xmin);
	double dy = 0.04*(ysize);
	double tmax_x[] = {x0, x0 + dx, x0 - dx, x0};
	double tmax_y[] = {y0, y0 + dy, y0 + dy, y0};
	tmax = new TPolyLine(4,tmax_x,tmax_y);
	tmax->SetLineColor(kRed);
	tmax->SetFillColor(kRed);
	tmax->Draw("f");

	// write mode location and shaded band
	// format of the number
	double order = TMath::Log10(TMath::Abs(x0));
	char sf='f';
	if ( order>6 || order<-3 )
		sf='e';

	TLatex * tmax_text = new TLatex();
	tmax_text->SetTextSize(0.045);
	tmax_text->SetTextFont(62);
	tmax_text->SetTextAlign(22); // center

	double xprint=(xmax+xmin)/2.;
	double yprint=ysize*(1-1.4*tmax_text->GetTextSize());

	if(TMath::Abs(limit)<50) // just to ensure there's some sense in the number
		tmax_text->DrawLatex(xprint,yprint,
			Form( Form("%%s* = %%%c ^{+%%%c}_{ -%%%c}",sf,sf,sf),
				fHistogram->GetXaxis()->GetTitle(), x0, max-x0, x0-min));
	else if (limit<0)
		tmax_text->DrawLatex(xprint,yprint,
			Form( Form("%%s* (%d%%%% CL) < %%%c",-limit,sf),
				fHistogram->GetXaxis()->GetTitle(), max));
	else if (limit>0)
		tmax_text->DrawLatex(xprint,yprint,
			Form( Form("%%s* (%d%%%% CL) > %%%c",limit,sf),
				fHistogram->GetXaxis()->GetTitle(), min));

}

// --------------------------------------------------------- 
 
void BCH1D::GetSmallestInterval(double & min, double & max, double content)
{
	int nbins = fHistogram -> GetNbinsX();
	int mins1 = -1;
	int mine1 = -1;
	double mininterval = nbins;

	// loop over start value
	for (int s1 = 1; s1 <= nbins; s1++)
	{
		// loop over end value
		for (int e1 = s1; e1 <= nbins; e1++)
		{
			double integral = fHistogram -> Integral(s1, e1, "width");
			double left  = 0.0;
			double right = 0.0;

			if (s1 > 1)
				left = fHistogram -> GetBinContent(s1 - 1) * fHistogram -> GetBinWidth(s1 - 1);
			if (e1 < nbins)
				right = fHistogram -> GetBinContent(e1 + 1) * fHistogram -> GetBinWidth(e1 + 1);

			// check to the right
			if ((integral < content) &&
				((integral + right) >= content) &&
				((e1 - s1 + 2) < mininterval))
			{
				mininterval = e1 - s1 + 2;
				mins1 = s1;
				mine1 = e1 + 1;
			}

			if (s1 == e1 && integral >= content)
			{
				mininterval = 1;
				mins1 = s1;
				mine1 = e1;
			}
		}
	}

	min = fHistogram->GetBinCenter(mins1);
	max = fHistogram->GetBinCenter(mine1);
}

// --------------------------------------------------------- 
