#include "BCLog.h"

#include "style.c"
#include "BCBenchmarkMCMC.h"
#include <TF1.h>

// ---------------------------------------------------------

int main()
{
	SetStyle();
	BCLog::OpenLog("log.txt", BCLog::detail, BCLog::detail);

	// function range
	double xmin =  0.0;
	double xmax = 50.0;

	// define test functions
	TF1 * testfunction2 = new TF1("testfunction2", "x*x*sin(x)*sin(x)", xmin, xmax);

	// number of runs with different number of iterations
	int n = 5;
	// factor to increase the number of iterations in a subsequent run
	int step = 10;

	// arrays for storing the chi2 and number of iterations for later plotting
	double * x = new double[n];
	double * y = new double[n];

	// do the runs
	int f=1;
	for(int i=0;i<n;i++)
	{
		// create model and set it up
		BCBenchmarkMCMC * benchmark = new BCBenchmarkMCMC("benchmark");
		benchmark -> AddParameter("x",  xmin, xmax);
		benchmark -> SetNbins(500, 0);
		benchmark -> SetTestFunction(testfunction2);

		f *= step;
		benchmark -> MCMCSetNIterationsRun(100*f);
		// run MCMC
		benchmark -> MarginalizeAll();

		// compare the generated distribution with the true function
		double chi2 = benchmark -> PerformTest(benchmark -> GetBestFitParameters(),
				0,
				benchmark -> GetMarginalized("x"),
				true,
				TString::Format("test-%d.ps",100*f));

		delete benchmark;

		x[i] = 100*f;
		y[i] = chi2;
	}

	// plot the Chi2 vs. number of iterations
	TCanvas * c = new TCanvas();
	gPad -> SetLogx();
	gPad -> SetLogy();
	TGraph * g_chi2 = new TGraph(n,x,y);
	g_chi2 -> SetMarkerStyle(20);
	g_chi2 -> SetMarkerSize(.5);
	g_chi2 -> Draw("apl");

	c -> Print("chi2.ps");


	delete[] x;
	delete[] y;
	delete g_chi2;
	delete c;

	return 0;
}

// ---------------------------------------------------------
