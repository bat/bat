#include "BCBenchmarkMCMC.h"
#include "BCMath.h"

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

	hist_prob -> Scale(normalization/hist_prob->Integral("width"));

	// initialize chi2
	double chi2 = 0;

	// get number of bins
	int nbins = hist_temp -> GetNbinsX();

	// loop over bins
	for (int i = 1; i <= nbins; ++i)
	{
		double y = hist_prob -> GetBinContent(i);
		double sigma2 = hist_temp -> GetBinContent(i);

		double min = hist_temp -> GetBinLowEdge(i);
		double max = min + hist_temp -> GetBinWidth(i);
		double y0 =  fTestFunction -> Integral(min,max);

		if (sigma2 > 0.)
			chi2 += (y-y0)*(y-y0)/2./sigma2;
	}

	// divide by the number of d.o.f.
	chi2 /= double(nbins);

	BCLog::Out(BCLog::summary,BCLog::summary,
			TString::Format(" Chi2 from test = %f (NDoF = %d)",chi2,nbins));
	BCLog::Out(BCLog::summary,BCLog::summary,
			TString::Format(" Results printed to file %s",filename));

	// print to file
	if (flag_print)
	{
		TCanvas * can = new TCanvas("can");
		can -> cd();
//		gPad -> SetLogy();
		hist_prob -> SetStats(0);
		hist_prob -> GetYaxis() -> SetTitle("f(x)");
		hist_prob -> Draw();
		fTestFunction -> SetNpx(1000);
		fTestFunction -> SetLineColor(kRed);
		fTestFunction -> SetLineWidth(1);
		fTestFunction -> Draw("c same");
		can -> Print(filename);
	}

	// free memory
	delete hist_prob;

	return chi2;
}

// ---------------------------------------------------------
