//
// This ROOT macro is part of BAT and can only be run if BAT
// was installed correctly. The macro shows an example of fitting
// an efficiency with a function defined by the user. The input are
// two histograms, one being a subset of the other. In the fit the
// uncertainties are considered to be poissonian.
//
// The macro can be run from within ROOT via commands
//
//    root[1] .x efficientyFitterExample.C
//
// or
//
//    root[1] .L efficientyFitterExample.C
//    root[2] efficientyFitterExample()
//
// or from the command line
//
//    $ root efficientyFitterExample.C
//


// The data fitted are generated randomly as a signal peak (gaussian)
// plus a flat background using a function CreateHistogram(nbins, ns, nb, seed)
// The arguments are 'nbins' - number of bins of the histogram,
// 'ns' - number of signal events, 'nb' - number of background events,
// 'seed' - initial seed for the random number generator. The location
// and the width of the signal peak can be set up using the variables
// 'mean' and 'sigma' below.

void CreateHistograms(int nbins, int nevents, TH1D * hist1, TH1D * hist2, int seed = 0);
double fitfunction(double * x, double * par);

const double mean  = 40.0;
const double sigma = 15.0;


// ---------------------------------------------------------
void efficientyFitterExample()
{
	int nbins = 100;

	// open log file
	BCLog::OpenLog("log.txt", BCLog::detail, BCLog::detail);

	// set style
	BCAux::SetStyle();

	// create new histograms
	hist1 = new TH1D("data1", ";x;N", nbins, 0.0, 100.0);
	hist1 -> SetStats(kFALSE);
	hist1 -> SetLineColor(kBlack);
	hist1 -> SetFillColor(5);
	hist1 -> SetFillStyle(1001);

	hist2 = new TH1D("data2", ";x;N", nbins, 0.0, 100.0);
	hist2 -> SetStats(kFALSE);
	hist2 -> SetLineColor(kRed);

	// create data
	CreateHistograms(nbins, 1000, hist1, hist2, 1000);

	// define a fit function
	TF1 * f1 = new TF1("f1", fitfunction, 0.0, 100.0, 2);
	f1 -> SetParLimits(0, 20.0, 70.0);
	f1 -> SetParLimits(1, 0.0,  20.0);

	// create a new efficiency fitter
	BCEfficiencyFitter * hef = new BCEfficiencyFitter(hist1, hist2, f1);

	// set options for evaluating the fit function
	hef -> SetFlagIntegration(false);

	// set options for MCMC
	hef -> MCMCSetNIterationsMax(100);
	hef -> MCMCSetNIterationsRun(100);

	// perform fit
	hef -> Fit();

	// print data and fit
	TCanvas * c = new TCanvas("c1","",800,400);
	c -> Divide(2, 1);
	c -> cd(1);
	hist1 -> Draw();
	hist2 -> Draw("same");
	c -> cd(2);
	hef -> DrawFit("", true); // draw with a legend

	// print to file
	c -> Print("fit.ps");

	// print marginalized distributions
	hef -> PrintAllMarginalized("distributions.ps");

}

// ---------------------------------------------------------
void CreateHistograms(int nbins, int nevents, TH1D * hist1, TH1D * hist2, int seed)
{
	// initialize random number generator
	gRandom = new TRandom3(seed);

	// fill histogram 1
	for (int i = 0; i < nevents ; ++i)
	{
		double x = gRandom -> Landau(30., 8.);
		hist1 -> Fill(x);
	}

	// fill histogram 2
	for (int i = 1; i <= nbins; ++i)
	{
		double x = hist1 -> GetXaxis() -> GetBinCenter(i);
		double eff = 0;
		if (x < mean)
			eff = TMath::Erfc((mean-x)/sigma)/2.;
		else
			eff = TMath::Erf((x-mean)/sigma)/2. + 0.5;
//		int n = gRandom -> Poisson(nevents);
//		hist1 -> SetBinContent(i, n);
//		hist2 -> SetBinContent(i, gRandom -> Binomial(n, eff));
		hist2 -> SetBinContent(i, gRandom -> Binomial(hist1 -> GetBinContent(i), eff));
	}
}

// ---------------------------------------------------------
double fitfunction(double * x, double * par)
{
	double ff;

	if (x[0] < par[0])
		ff = TMath::Erfc((par[0]-x[0])/par[1])/2.;
	else
		ff = TMath::Erf((x[0]-par[0])/par[1])/2. + 0.5;

	return ff;
}

// ---------------------------------------------------------
