
TGraphErrors * CreateGraph(int n, int seed = 0);

const double slope  = 1.0;
const double offset = 11.0;
const double sigma  = 5.0;

// ---------------------------------------------------------

void runBATexample()
{
	BCAux::SetStyle();

	// create data
	TGraphErrors * graph = CreateGraph(10, 1000);

	// define a fit function
	TF1 * f1 = new TF1("f1", "[0] + [1]*x", 0.0, 100.0);
	f1 -> SetParLimits(0, -20.0, 30.0);
	f1 -> SetParLimits(1,   0.5,  1.5);

	// create a new graph fitter
	BCGraphFitter * gf = new BCGraphFitter(graph, f1);

	// set options for MCMC
	gf -> MCMCSetNIterationsRun(10000);

	// perform fit
	gf -> Fit();

	// create histogram
	TH2D * hist_axes = new TH2D("hist_axes", ";x;y", 1, 0.0, 99.5, 1, 0.0, 120.0);
	hist_axes -> SetStats(kFALSE);

	// print data and fit
	TCanvas * c1 = new TCanvas("c1");
	c1 -> cd();
	hist_axes -> Draw();
	gf -> DrawFit("SAME", true); // draw on top and add a legend

	// print to file
	c1 -> Print("fit.ps");

	// print marginalized distributions
	gf -> PrintAllMarginalized("distributions.ps");
}

// ---------------------------------------------------------

TGraphErrors * CreateGraph(int n, int seed)
{
	// initialize random number generator
	gRandom = new TRandom3(seed);
	
	// define arrays
	double * x  = new double[n];
	double * y  = new double[n];
	double * ey = new double[n];

	// define x and y-values and the uncertainties on y
	for (int i = 0; i < n; ++i)
	{
		x[i] = 100. / double(n) * double(i+0.5);
		y[i] = gRandom -> Gaus(offset + slope*x[i], sigma);
		ey[i] = sigma;
	}

	// create new graph
	TGraphErrors * graph = new TGraphErrors(n, x, y, 0, ey);
	graph -> SetMarkerStyle(20);
	graph -> SetMarkerSize(1.5);

	// return the graph
	return graph;
}

// ---------------------------------------------------------
