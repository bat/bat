#include "BCH1D.h"
#include "BCMath.h" 

#include <TH2.h>
#include <TLine.h>
#include <TPolyLine.h>
#include <TLatex.h>
#include <TError.h>
#include <TCanvas.h>
#include <TMarker.h>
#include <TLegend.h>

// ---------------------------------------------------------

BCH1D::BCH1D()
{
	
	fHistogram = 0;
	fDefaultCLLimit = 95.; // in percent 

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

	// use ROOT function to calculat quantile. 

	fHistogram -> GetQuantiles(nquantiles, quantiles, probsum); 

	return quantiles[0]; 

}

// --------------------------------------------------------- 

double BCH1D::GetIntegral(double valuemin, double valuemax) 
{

	double integral = 0; 

	int binmin = fHistogram -> FindBin(valuemin); 
	int binmax = fHistogram -> FindBin(valuemax); 

	// use ROOT function to calculate integral. 

	integral = fHistogram -> Integral(binmin, binmax); 

	return integral; 

}

// --------------------------------------------------------- 

double BCH1D::GetPValue(double probability)
{

	// get maxmimum value of the distribution. 

	double valuemax = fHistogram -> GetBinCenter(fHistogram -> GetNbinsX()); 

	// use ROOT function to calculate the integral from "probability" 
	// to the maximum. 

	double integral = this -> GetIntegral(probability, valuemax); 

	return integral; 

}

// --------------------------------------------------------- 

void BCH1D::SetDefaultCLLimit(double limit)
{

	// check if limit is between 68% and 100%. Give a warning if not ... 

	if(limit>=100. || limit<68.)
		BCLog::Out(BCLog::warning,BCLog::warning,
			Form("BCH1D::SetDefaultCLLimit. Value %.1f out of range. Keeping default %.1f%% CL limit.",limit,fDefaultCLLimit));

	// ... or set value if true. 

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

	// define temporary line. 

	TLine * line;

	// check drawing option. 

	switch(options)
		{
			// Draw a band between 16% and 84% probability. 
			// If the mode is outside the band only draw a limit. 

		case 0:

			if (fabs(ovalue) >= 100 || fabs(ovalue) < 68)
				{
					min = this -> GetQuantile(.16);
					max = this -> GetQuantile(.84);

					if ( mode > max || fHistogram -> FindBin(mode) == fHistogram -> GetNbinsX() )
						{
							min = this -> GetQuantile(1.-(double)fDefaultCLLimit/100.);
							max = fHistogram->GetXaxis()->GetXmax();
							ovalue = fDefaultCLLimit;
						}
					else if (mode < min || fHistogram->FindBin(mode)==1)
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

			// Draw a line at "ovalue". 

		case 1:

			fHistogram ->Draw();
			line = new TLine();
			line -> SetLineColor(kRed);
			line -> DrawLine(ovalue, 0., ovalue, 1.05 * fHistogram -> GetMaximum());
		
			break;

			// Draw a shaded band at the smallest interval. 

		case 2:

			if(ovalue<50) // just to ensure there's some sense in the number
				ovalue = 68.; // default is 68%

			this->GetSmallestInterval(min, max, ovalue/100.);
			this->PrintShadedLimits(mode, min, max, 0.);

			break;

			// Sort out bad options and warn. 

		default:

			BCLog::Out(BCLog::warning, BCLog::warning, Form("BCH1D::Print. Invalid option %d",options));
			break;
		}

	// print to file. 

	canvas -> Print(filename);
}

// --------------------------------------------------------- 

void BCH1D::PrintShadedLimits(double mode, double min, double max, double limit)
{

	double maximum = fHistogram -> GetMaximum();

	double x0 = mode;
	double y0 = fHistogram->GetBinContent( fHistogram->FindBin(mode) );

	double ysize = maximum*1.2;

	double xmin = fHistogram->GetXaxis()->GetXmin();
	double xmax = fHistogram->GetXaxis()->GetXmax();

	// draw histogram with axes first

	TH2D * hsc = new TH2D(Form("h2scale_%s_%d",fHistogram->GetName(),BCLog::GetHIndex()),"",
												10, xmin, xmax, 10, 0., ysize);
	hsc -> SetStats(0);
	hsc -> SetXTitle(fHistogram->GetXaxis()->GetTitle());
	hsc -> SetYTitle(fHistogram->GetYaxis()->GetTitle());
	hsc -> Draw();

	// draw histogram
	
	fHistogram -> Draw("same");

	// draw yellow shaded region between min and max

	TH1D * hist_shaded = this->GetSubHisto(min,max,Form("%s_sub_%d",fHistogram->GetName(),BCLog::GetHIndex()));
	hist_shaded -> SetFillStyle(1001);
	hist_shaded -> SetFillColor(kYellow);

	// draw shaded histogram
	
	hist_shaded -> Draw("same");

	gPad->RedrawAxis();

	// draw line and triangle for mode

	TLine * line;
	TPolyLine * tmax;

	if(fabs(limit)<50) // just to ensure there's some sense in the number
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

	double order = log10(fabs(x0));
	char sf='f';
	if ( order>6 || order<-3 )
		sf='e';

	TLatex * tmax_text = new TLatex();
	tmax_text->SetTextSize(0.045);
	tmax_text->SetTextFont(62);
	tmax_text->SetTextAlign(22); // center

	double xprint=(xmax+xmin)/2.;
	double yprint=ysize*(1-1.4*tmax_text->GetTextSize());

	if(fabs(limit)<50) // just to ensure there's some sense in the number
		tmax_text->DrawLatex(xprint,yprint,
			Form( Form("%%s* = %%%c ^{+%%%c}_{ -%%%c}",sf,sf,sf),
				fHistogram->GetXaxis()->GetTitle(), x0, max-x0, x0-min));

	else if (limit<0)
		tmax_text->DrawLatex(xprint,yprint,
			Form( Form("%%s* (%d%%%% CL) < %%%c",-(int)limit,sf),
				fHistogram->GetXaxis()->GetTitle(), max));

	else if (limit>0)
		tmax_text->DrawLatex(xprint,yprint,
			Form( Form("%%s* (%d%%%% CL) > %%%c",(int)limit,sf),
				fHistogram->GetXaxis()->GetTitle(), min));

}

// ---------------------------------------------------------

void BCH1D::GetSmallestInterval(double & min, double & max, double content)
{

	if(content<=0. || content>= 1.)
		return;

	int nbins = fHistogram -> GetNbinsX();

	double factor = fHistogram->Integral("width");
	if(factor==0)
		return;

	// normalize if not yet done

	fHistogram->Scale(1./factor);

	double xmin = fHistogram->GetXaxis()->GetXmin();
	double xmax = fHistogram->GetXaxis()->GetXmax();
	double width = xmax - xmin;

	double xdown=xmin;
	double xup=xmax;

	int ndiv = 100;
	if(nbins<100)
		ndiv = 1000;
	if(nbins>1000)
		ndiv = 10;

	int warn=0;

	// loop through the bins

	for(int i=1;i<nbins+1;i++)
	{
		if(fHistogram->Integral(i,nbins,"width") < content)
			break;

		// get width of the starting bin

		double firstbinwidth = fHistogram->GetBinWidth(i);

		// divide i-th bin in ndiv divisions
		// integrate starting at each of these divisions

		for(int j=0;j<ndiv;j++)
		{
			double dxdown = (double)(ndiv-j)/(double)ndiv * firstbinwidth;
			xdown = fHistogram->GetBinLowEdge(i) + firstbinwidth - dxdown;
			double integral = dxdown * fHistogram->GetBinContent(i);

			if(integral>content)
			{
				// content fits within one bin

				xup = xdown + content / fHistogram->GetBinContent(i);
				warn = 1;
			}
			else
			{
				for(int k=i+1;k<nbins+1;k++)
				{
					double thisbin = fHistogram->GetBinContent(k) * fHistogram->GetBinWidth(k);

					if(integral+thisbin==content)
					{
						xup = fHistogram->GetBinLowEdge(k+1);
						break;
					}

					if(integral+thisbin>content)
					{
						xup = fHistogram->GetBinLowEdge(k) + (content-integral)/thisbin * fHistogram->GetBinWidth(k);
						integral += thisbin * (xup-fHistogram->GetBinLowEdge(k))/fHistogram->GetBinWidth(k);
						break;
					}

					integral += thisbin;
				}
			}
			if(integral < content)
				continue;
			if(xup - xdown < width)
			{
				// store necessary information
				width = xup - xdown;
				xmin  = xdown;
				xmax  = xup;
			}
		}
	}

	if(warn)
	{
		BCLog::Out(BCLog::warning,BCLog::warning,
			Form("BCH1D::GetSmallestInterval. The requested content of %d%% fits within one bin.",(int)(content*100)));
		BCLog::Out(BCLog::warning,BCLog::warning,
			"BCH1D::GetSmallestInterval. MAKE FINER BINNING (OR CHANGE RANGE) !!!");
	}

	// restore normalization to state before calling this method

	fHistogram->Scale(factor);

	min=xmin;
	max=xmax;

}

// ---------------------------------------------------------

double BCH1D::IntegralWidth(double min, double max)
{
	int imin = fHistogram->FindBin(min);
	int imax = fHistogram->FindBin(max);

	int nbins = fHistogram->GetNbinsX();

	// if outside of histogram range, return -1.
	if ( imin<1 || imin>nbins || imax<1 || imax>nbins )
		return -1.;

	if ( imin==imax )
		return -1.;

	// swap values if necessary
	if (imin>imax)
	{
		int i=imin;
		double x=min;
		imin=imax, min=max;
		imax=i, max=x;
	}

	// calculate first bin
	double first = ( fHistogram->GetBinLowEdge(imin+1) - min ) * fHistogram->GetBinContent(imin);

	// calculate last bin
	double last = ( max - fHistogram->GetBinLowEdge(imax) ) * fHistogram->GetBinContent(imax);

	// calculate inbetween
	double inbetween=0.;
	if(imax-imin>1)
		inbetween = fHistogram->Integral(imin+1, imax-1, "width");

	return first + inbetween + last;
}

// ---------------------------------------------------------

TH1D * BCH1D::GetSubHisto(double min, double max, const char * name)
{
	if(min>max)
	{
		double t=min;
		min=max;
		max=t;
	}

	int imin = fHistogram->FindBin(min);
	int imax = fHistogram->FindBin(max);

	int nbins = fHistogram->GetNbinsX();
	double xmin = fHistogram->GetXaxis()->GetXmin();
	double xmax = fHistogram->GetXaxis()->GetXmax();

	if( min==max || (min<=xmin && max>=xmax) )
	{
		TH1D * h0 = (TH1D*) fHistogram->Clone(name);
		return h0;
	}

	if (imin<1)
	{
		imin=1;
		min=xmin;
	}
	if (imax>nbins)
	{
		imax=nbins;
		max=xmax;
	}

	double * xb = new double[nbins+3]; // nbins+1 original bin edges + 2 new bins
	int n=0; // counter

	int domin=1;
	int domax=1;

	for (int i=1;i<nbins+2;i++)
	{
		double x0 = fHistogram->GetBinLowEdge(i);

		if (min<x0 && domin)
		{
			xb[n++]=min;
			domin=0;
		}
		else if (min==x0)
			domin=0;

		if (max<x0 && domax)
		{
			xb[n++]=max;
			domax=0;
		}
		else if (max==x0)
			domax=0;

		xb[n++]=x0;

	}

	// now define the new histogram
	TH1D * h0 = new TH1D(name,"",n-1,xb);
	for(int i=1;i<n;i++)
	{
		double x0 = h0->GetBinCenter(i);
		if(x0<min || x0>max)
			continue;

		int bin=fHistogram->FindBin(x0);
		h0->SetBinContent(i, fHistogram->GetBinContent(bin));
	}

	return h0;
}

// ---------------------------------------------------------
