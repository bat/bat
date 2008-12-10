#include <BCHistogramFitter.h>
#include <BCLog.h>
#include <BCModelOutput.h>

#include <TH1D.h>
#include <TF1.h>

#include <TMath.h>

#include "style.c"

// ---------------------------------------------------------

int main()
{

	// ---------------------------------------------------------
	// set style
	// ----------------------------------------------------------
	SetStyle();

	// ---------------------------------------------------------
	// open log file
	// ---------------------------------------------------------
	BCLog::OpenLog("log.txt", BCLog::detail, BCLog::detail);

	// ---------------------------------------------------------
	// create a histogram with some data
	// ---------------------------------------------------------
	TH1D * hist = new TH1D("histogram", ";E [GeV];N", 20, 0., 100.);
	hist -> SetStats(kFALSE);

	gRandom = new TRandom(1000);

	// fill "signal" events
	for (int i = 0; i < 50; ++i)
		hist -> Fill(gRandom -> Gaus(50., 10.));

	// fill "background" events
	for (int i = 0; i < 100; ++i)
		hist -> Fill(100.0 * gRandom -> Uniform());

	// ---------------------------------------------------------
	// define a fit function
	// ---------------------------------------------------------
	TF1 * func = new TF1("func", "[0] / sqrt(2.0 * 3.1416) / [2] * exp(-(x-[1])*(x-[1])/2./[2]/[2]) + [3]", 0., 100.);
	func -> SetParLimits(0, 0.0, 200.);
	func -> SetParLimits(1, 0.0, 100.);
	func -> SetParLimits(2, 0.1, 25.0);
	func -> SetParLimits(3, 0.0, 2.0);

	// ---------------------------------------------------------
	// initialize histogram fitter and perform the fit 
	// ---------------------------------------------------------
	BCHistogramFitter * hf = new BCHistogramFitter(); 
	hf -> Fit(hist, func); 

	hf -> PrintAllMarginalized("plots.ps");

	// ---------------------------------------------------------
	// print histogram
	// ---------------------------------------------------------
	TCanvas * c1 = new TCanvas();
//	c1 -> cd();
	hist -> Draw();
	hist -> GetYaxis() -> SetRangeUser(0., 22.);
	hf -> GetErrorBandGraph(0.16, 0.84) -> Draw("f same");
	hist -> Draw("SAME");

	// go from a number density to a number of expected events
	func -> SetParameter(0, func -> GetParameter(0) * hist -> GetBinWidth(0));
	func -> SetParameter(3, func -> GetParameter(3) * hist -> GetBinWidth(0));
	func -> Draw("SAME");

	c1 -> Print("data.ps");

	// ---------------------------------------------------------
	// close log file
	// ---------------------------------------------------------
	BCLog::CloseLog();

	return 0;

}

// ---------------------------------------------------------

