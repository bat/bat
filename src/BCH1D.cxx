#include "BCH1D.h"

#include <TCanvas.h>
#include <TLine.h>
#include <TMarker.h>
#include <TLegend.h>

// ---------------------------------------------------------

BCH1D::BCH1D()
{
	fHistogram = 0;

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
  
  double valuemax = fHistogram -> GetBinCenter(fHistogram -> GetNbinsX()); 
  
  double integral = this -> GetIntegral(probability, valuemax); 

  return integral; 

}

// --------------------------------------------------------- 

void BCH1D::Print(char* filename, int options, double value) 
{
	// create temporary canvas
	TCanvas * canvas = new TCanvas();
	canvas -> cd();

	double maximum = fHistogram -> GetMaximum();

	double location16 = this -> GetQuantile(.16);
	double location84 = this -> GetQuantile(.84);

	switch(options)
	{
		// draw lines
		case 1:
		{
			// draw histogram
			fHistogram -> Draw();

			TLine * line = new TLine();
			// draw line for mode
			line -> DrawLine(this -> GetMode(), 0., this -> GetMode(), 1.05 * maximum);
			// draw line for 16%
			line -> DrawLine(location16, 0., location16, 1.05 * maximum);
			// draw line for 84%
			line -> DrawLine(location84, 0., location84, 1.05 * maximum);

			break;
		}

		// draw line at value
		case 2:
		{
			// draw histogram
			fHistogram -> Draw();

			TLine * line = new TLine();
			line -> SetLineColor(kRed);
			line -> DrawLine(value, 0., value, 1.05 * maximum);
			break;
		}

		// nice printing
		// draw 68% band and print mode
		case 3:
		{
			double x0 = this->GetMode();
			double y0 = maximum;

			double ysize = maximum*1.2;

			// draw yellow shaded region for the 16% and 84% quantiles with a red triangle on top of the mode
			TH1D * hist_shaded =  (TH1D *) fHistogram -> Clone();
			hist_shaded -> SetFillStyle(1001);
			hist_shaded -> SetFillColor(kYellow);

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

			// remove entries on the sides
			int nbins = hist_shaded -> GetNbinsX();
			int lowerborder = hist_shaded -> FindBin(this -> GetQuantile(.16));
			int upperborder = TMath::Min(hist_shaded -> FindBin(this -> GetQuantile(.84)) + 1, nbins);

			for (int i = 1; i < lowerborder; i++)
				hist_shaded -> SetBinContent(i, 0.);

			for (int i = upperborder; i <= nbins; i++)
				hist_shaded -> SetBinContent(i, 0.);

			// draw shaded histogram
			hist_shaded -> Draw("same");

			gPad->RedrawAxis();

			// draw line for mode
			TLine * line = new TLine();
			line -> SetLineStyle(2);
			line -> DrawLine(this -> GetMode(), 0., this -> GetMode(), maximum);

			// draw triangle
			double dx = 0.01*(xmax-xmin);
			double dy = 0.04*(ysize);
			double tmax_x[] = {x0, x0 + dx, x0 - dx, x0};
			double tmax_y[] = {y0, y0 + dy, y0 + dy, y0};
			TPolyLine * tmax = new TPolyLine(4,tmax_x,tmax_y);
			tmax->SetLineColor(kRed);
			tmax->SetFillColor(kRed);
			tmax->Draw("f");

			// write mode location and 68% band
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
			tmax_text->DrawLatex(xprint,yprint,
				Form( Form("%%s* = %%%c ^{+%%%c}_{ -%%%c}",sf,sf,sf),
					fHistogram->GetXaxis()->GetTitle(), x0, location84-x0, x0-location16));

	      break;
		}

		// bad options
		default:
			BCLog::Out(BCLog::warning, BCLog::warning, Form("BCH1D::Print. Invalid option %d",options));
			break;
	}

	// print to file 
	canvas -> Print(filename); 
}

// --------------------------------------------------------- 

void BCH1D::Print(char* filename, int options) 
{
  
  double value = 0; 

  this -> Print(filename, options, value); 

}

// --------------------------------------------------------- 
 
void BCH1D::Print(char* filename) 
{
  
  int options = 0; 
  double value = 0; 

  this -> Print(filename, options, value); 

}

// --------------------------------------------------------- 
 





