//
// This ROOT macro is part of BAT and can only be run if BAT
// was installed correctly. The macro shows an example of fitting
// a histogram using a function defined by the user. Poissonian
// uncertainties are used in the fit.
//
// The macro can be run from within ROOT via commands
//
//    root[1] .x hFitterExample.C
//
// or
//
//    root[1] .L hFitterExample.C
//    root[2] hFitterExample()
//
// or from the command line
//
//    $ root hFitterExample.C
//

#include <TH1D.h>
#include <TF1.h>
#include <TRandom.h>
#include <TCanvas.h>

void hFitterExample()
{
	// set up nice style for plots
	BCAux::SetStyle();

	// uncomment the this command if you want more printouts during the run
//	BCLog::SetLogLevel(BCLog::detail);

	// create a histogram with some data
	TH1D * hist = new TH1D("histogram", ";E [GeV];N", 20, 0., 100.);
	hist -> SetStats(kFALSE);

	gRandom = new TRandom(1000);

	// fill "signal" events
	for (int i = 0; i < 50; ++i)
		hist -> Fill(gRandom -> Gaus(50., 10.));

	// fill "background" events
	for (int i = 0; i < 100; ++i)
		hist -> Fill(100.0 * gRandom -> Uniform());

	// define a fit function
	TF1 * func = new TF1("func", "[0] / sqrt(2.0 * 3.1416) / [2] * exp(-(x-[1])*(x-[1])/2./[2]/[2]) + [3]", 0., 100.);
	func -> SetParLimits(0, 0.,  200.);
	func -> SetParLimits(1, 0.,  100.);
	func -> SetParLimits(2, 0.1,  25.);
	func -> SetParLimits(3, 0.,    2.);

	// initialize histogram fitter and perform the fit
	BCHistogramFitter * hf = new BCHistogramFitter();
	hf -> Fit(hist, func);

	// draw fit including the error band
	TCanvas * c1 = new TCanvas();
	hf -> DrawFit();
//	c1 -> Print("data.ps");

	// uncomment the next command if you want to
	// draw marginalized posterior distributions for all parameters
	// and print them into a postscript file
//	hf -> PrintAllMarginalized("plots.ps");
}
