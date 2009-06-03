#include "BAT/BCMath.h"
#include "BAT/BCH2D.h"
#include "BAT/BCLog.h"

#include <TROOT.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TH2D.h>

#include "BCBenchmarkMCMC2D.h"

// ---------------------------------------------------------

BCBenchmarkMCMC2D::BCBenchmarkMCMC2D(const char* name) : BCModel(name)
{}

// ---------------------------------------------------------

double BCBenchmarkMCMC2D::LogLikelihood(std::vector <double> parameters)
{return log (fTestFunction -> Eval(parameters[0],parameters[1]));}

// ---------------------------------------------------------

double BCBenchmarkMCMC2D::PerformTest(
		std::vector<double> parameters,
		int index,
		BCH2D * hist,
		bool flag_print,
		const char * filename)
{
	// get histogram from BCH2D and clone it
	TH2D * hist_temp = hist -> GetHistogram();
	TH2D * hist_prob = (TH2D*) hist_temp -> Clone();

	double xmin = this->GetParameter(0)->GetLowerLimit();
	double xmax = this->GetParameter(0)->GetUpperLimit();
	double ymin = this->GetParameter(1)->GetLowerLimit();
	double ymax = this->GetParameter(1)->GetUpperLimit();
	double normalization = fTestFunction -> Integral(xmin,xmax,ymin,ymax);

	hist_prob -> Sumw2();
	hist_prob -> Scale(normalization/hist_prob->Integral("width"));

	// initialize chi2
	double chi2 = 0;

	// get number of bins
	int nbinsx = hist_prob -> GetNbinsX();
	int nbinsy = hist_prob -> GetNbinsY();
	double XbinWidth = (xmax-xmin)/nbinsx;
	double YbinWidth = (ymax-ymin)/nbinsy;

	// loop over bins
	TH2D* hdiff = new TH2D("hdiff","",nbinsx,xmin,xmax,nbinsy,ymin,ymax);
	for (int i = 1; i <= nbinsx; ++i) {
		for (int j = 1; j <= nbinsy; ++j) {
			double v = hist_prob -> GetBinContent(i,j);
			double fxmin = xmin+(i-1)*XbinWidth;
			double fxmax = fxmin+XbinWidth;
			double fymin = ymin+(j-1)*YbinWidth;
			double fymax = fymin+YbinWidth;
			double v0 =  fTestFunction -> 
				Integral(fxmin,fxmax,fymin,fymax)/(XbinWidth*YbinWidth);
			double vmv0 = v - v0;

			double dv = hist_prob -> GetBinError(i,j);
			if (dv > 0.) {
				chi2 += (vmv0*vmv0)/(dv*dv);
				hdiff->SetBinContent(i,j,vmv0/dv);
			}
		}
	}

	// divide by the number of d.o.f.
	chi2 /= double(nbinsx*nbinsy);

	BCLog::Out(BCLog::summary,BCLog::summary,
			TString::Format(" Chi2 from test = %f (NDoF = %d)",chi2,nbinsx*nbinsy));

	// print to file
	if (flag_print)
	{
		TCanvas *can = (TCanvas*) gROOT->GetListOfCanvases()->FindObject("can");
		if (!can)  can = new TCanvas("can");
		can -> Clear();
		//gPad -> SetLogz();
		gStyle -> SetOptStat(0);
		hdiff -> GetXaxis() -> SetTitle("x");
		hdiff -> GetYaxis() -> SetTitle("y");
		hdiff -> Draw("colz");
		TLatex tx2;
		tx2.SetNDC();
		tx2.DrawLatex(0.3,0.9,Form("#chi^{2}/NDF: %f",chi2));
		can -> Print(filename);
	}
	BCLog::Out(BCLog::summary,BCLog::summary,
			TString::Format(" Results printed to file %s",filename));

	// free memory
	delete hist_prob;
	delete hdiff;

	return chi2;
}

// ---------------------------------------------------------
