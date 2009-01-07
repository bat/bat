void CreateHistograms(int nbins, int nevents, TH1D * hist1, TH1D * hist2, int seed = 0); 
double fitfunction(double *x, double *par); 

const double mean  = 40.0; 
const double sigma = 15.0;

// ---------------------------------------------------------
#include "../../include/style.c"
// ---------------------------------------------------------

int runBATexample()
{
	int nbins = 100;

	// open log file
	BCLog::OpenLog("log.txt", BCLog::detail, BCLog::detail);

	// set style
	SetStyle(); 

	// create new histograms 
	hist1 = new TH1D("data1", ";x [au];N", nbins, 0.0, 100.0); 
	hist1 -> SetStats(kFALSE); 
	hist1 -> SetLineColor(kBlack); 

	hist2 = new TH1D("data2", ";x [au];N", nbins, 0.0, 100.0); 
	hist2 -> SetStats(kFALSE); 
	hist2 -> SetLineColor(kRed); 

	// create data 
	CreateHistograms(nbins, 500, hist1, hist2, 1000); 

	// define a fit function 
	TF1 * f1 = new TF1("f1", fitfunction, 0.0, 100.0, 2); 
	f1 -> SetParLimits(0, 20.0, 70.0);
	f1 -> SetParLimits(1, 0.0,  20.0); 

	// create a new histogram fitter
	BCEfficiencyFitter * hef = new BCEfficiencyFitter(hist1, hist2, f1); 

	// set options for evaluating the fit function
	hef -> SetFlagIntegration(false); 

	// set options for MCMC 
	hef -> MCMCSetNIterationsMax(10000);
	hef -> MCMCSetNIterationsRun(10000);

	// perform fit 
	hef -> Fit();

	// print data and fit 
	TCanvas * c1 = new TCanvas("c1"); 
	c1 -> Divide(2, 1); 
	c1 -> cd(1);
	hist1 -> Draw(); 
	hist2 -> Draw("SAME"); 
	c1 -> cd(2); 
	hef -> DrawFit("", true); // draw with a legend

	// print to file 
	c1 -> Print("fit.ps"); 

	// print marginalized distributions 
	hef -> PrintAllMarginalized("distributions.ps"); 

	return 1;
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
	
	// fill histograms 
	for (int i = 1; i <= nbins; ++i)
		{
			double x = hist1 -> GetXaxis() -> GetBinCenter(i); 
			double eff = 0; 
			if (x < mean)
				eff = TMath::Erfc((mean-x)/sigma)/2.; 
			else
				eff = TMath::Erf((x-mean)/sigma)/2. + 0.5; 
			//			int n = gRandom -> Poisson(nevents); 
			//			hist1 -> SetBinContent(i, n);
			//			hist2 -> SetBinContent(i, gRandom -> Binomial(n, eff));
			hist2 -> SetBinContent(i, gRandom -> Binomial(hist1 -> GetBinContent(i), eff));
		}
}

// ---------------------------------------------------------

double fitfunction(double *x, double *par)
{
	double ff = 0;

	if (x[0] < par[0])
		ff = TMath::Erfc((par[0]-x[0])/par[1])/2.; 
	else
		ff = TMath::Erf((x[0]-par[0])/par[1])/2. + 0.5; 

	return ff; 

}

// ---------------------------------------------------------
