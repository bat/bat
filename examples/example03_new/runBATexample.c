TH1D * CreateHistogram(int nbins, int ns, int nb, int seed = 0);

const double mean  = 65.0; 
const double sigma =  5.0;

// ---------------------------------------------------------

int runBATexample()
{
	// create data 
	TH1D * hist = CreateHistogram(20, 100, 100); 

	// define a fit function 
	TF1 * f1 = new TF1("f1", "[0] / sqrt(2.0 * 3.1416) / [2] * exp(-(x-[1])*(x-[1])/2./[2]/[2]) + [3]", 0.0, 100.0); 
	f1 -> SetParLimits(0,  0.0, 200.0);
	f1 -> SetParLimits(1, 50.0,  90.0); 
	f1 -> SetParLimits(2,  0.1,  10.0); 
	f1 -> SetParLimits(3,  0.0,   2.0); 

	// create a new histogram fitter
	BCHistogramFitter * hf = new BCHistogramFitter(hist, f1); 

	// set options for MCMC 
	hf -> MCMCSetNIterationsRun(10000);

	// perform fit 
	hf -> Fit();

	// print data and fit 
	TCanvas * c1 = new TCanvas("c1"); 
	c1 -> cd(); 	
	hf -> DrawFit("", true); // draw with a legend

	// print to file 
	c1 -> Print("fit.ps"); 

	// print marginalized distributions 
	hf -> PrintAllMarginalized("distributions.ps"); 

	return 1;
}

// ---------------------------------------------------------

TH1D * CreateHistogram(int nbins, int ns, int nb, int seed)
{
	// initialize random number generator 
	gRandom = new TRandom3(seed); 

	// create new histogram 
	TH1D * hist = new TH1D("data", ";x [au];N", nbins, 0.0, 100.0); 
	hist -> SetStats(kFALSE); 
	
	// fill signal 
	for (int i = 0; i < ns; ++i)
		hist -> Fill(gRandom -> Gaus(mean, sigma)); 

	// fill background 
	for (int i = 0; i < nb; ++i)
		hist -> Fill(100.0 * gRandom  -> Uniform()); 

	// return the histogram
	return hist; 
}

// ---------------------------------------------------------
