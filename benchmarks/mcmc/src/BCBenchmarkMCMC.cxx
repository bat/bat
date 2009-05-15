#include "BAT/BCMath.h"
#include "BAT/BCH1D.h"
#include "BAT/BCLog.h"

#include <TROOT.h>
#include <TCanvas.h>
#include <TLatex.h>

#include "BCBenchmarkMCMC.h"

// ---------------------------------------------------------

BCBenchmarkMCMC::BCBenchmarkMCMC(const char* name) : BCModel(name)
{}

// ---------------------------------------------------------

double BCBenchmarkMCMC::LogLikelihood(std::vector <double> parameters)
{
	return log (fTestFunction -> Eval(parameters[0]));
}

// ---------------------------------------------------------

double BCBenchmarkMCMC::PerformTest(
		std::vector<double> parameters,
		int index,
		BCH1D * hist,
		bool flag_print,
		const char * filename)
{
	// get histogram from BCH1D and clone it
	TH1D * hist_temp = hist -> GetHistogram();
	TH1D * hist_prob = (TH1D*) hist_temp -> Clone();

	double xmin = this->GetParameter(0)->GetLowerLimit();
	double xmax = this->GetParameter(0)->GetUpperLimit();
	double normalization = fTestFunction -> Integral(xmin,xmax);

	hist_prob -> Sumw2();
	hist_prob -> Scale(normalization/hist_prob->Integral("width"));

	// initialize chi2
	double chi2 = 0;

	// get number of bins
	int nbins = hist_temp -> GetNbinsX();

	// loop over bins
	for (int i = 1; i <= nbins; ++i)
	{
		double y = hist_prob -> GetBinContent(i);
		double dy = hist_prob -> GetBinError(i);

		double min = hist_prob -> GetBinLowEdge(i);
		double max = min + hist_prob -> GetBinWidth(i);
		double y0 =  fTestFunction -> Integral(min,max)/(max-min);

		if (dy > 0.)
			chi2 += (y-y0)*(y-y0)/(dy*dy);
	}

	// divide by the number of d.o.f.
	chi2 /= double(nbins);

	BCLog::Out(BCLog::summary,BCLog::summary,
			TString::Format(" Chi2 from test = %f (NDoF = %d)",chi2,nbins));

	// print to file
	if (flag_print)
	{
		TCanvas *can = (TCanvas*) gROOT->GetListOfCanvases()->FindObject("can");
		if (!can)  can = new TCanvas("can");
		can -> Clear();
//		gPad -> SetLogy();
		hist_prob -> SetStats(0);
		hist_prob -> GetYaxis() -> SetTitle(fTestFunction->GetName());
		hist_prob -> Draw();
		fTestFunction -> SetNpx(1000);
		fTestFunction -> SetLineColor(kRed);
		fTestFunction -> SetLineWidth(2);
		fTestFunction -> Draw("c same");
		TLatex tx2;
		tx2.SetNDC();
		tx2.DrawLatex(0.3,0.9,Form("#chi^{2}/NDF: %f",chi2));
		can -> Print(filename);
	}
	BCLog::Out(BCLog::summary,BCLog::summary,
			TString::Format(" Results printed to file %s",filename));

	// free memory
	delete hist_prob;

	return chi2;
}

// ---------------------------------------------------------
