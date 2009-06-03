// test MCMC sampling for 2D functions
// Author: Jing Liu
//
#include "BAT/BCAux.h"
#include "BAT/BCLog.h"

#include "BCBenchmarkMCMC2D.h"

#include <TF2.h>

// ---------------------------------------------------------

int main()
{
	BCAux::SetStyle();

	BCLog::OpenLog("SamplingTest2D.log", BCLog::detail, BCLog::detail);

	BCBenchmarkMCMC2D* benchmark = new BCBenchmarkMCMC2D("2D Sampling Test");

	// add parameters
	double xmin =  4.0;
	double xmax = 16.0;
	double ymin =  4.0;
	double ymax = 16.0;
	benchmark -> AddParameter("x", xmin, xmax);
	benchmark -> AddParameter("y", xmin, xmax);




	// define and set test functions: flat plate
	TF2 * flatPlate = new TF2("flatPlate","10*x/x*y/y",xmin,xmax,ymin,ymax);
	benchmark -> SetTestFunction(flatPlate);

	// perform marginalization
	benchmark -> MCMCSetNIterationsRun(1000000);
	benchmark -> MarginalizeAll();

	// perform test
	benchmark -> PerformTest(
			benchmark -> GetBestFitParameters(),
			0,
			benchmark -> GetMarginalized("x","y"),
			true,
			"flatPlate.eps");




	// define and set test functions: x*y
	TF2 * xxy = new TF2("xxy","x*y",xmin,xmax,ymin,ymax);
	benchmark -> SetTestFunction(xxy);

	// perform marginalization
	benchmark -> MCMCSetNIterationsRun(1000000);
	benchmark -> MarginalizeAll();

	// perform test
	benchmark -> PerformTest(
			benchmark -> GetBestFitParameters(),
			0,
			benchmark -> GetMarginalized("x","y"),
			true,
			"xxy.eps");




	// define and set test functions: sin(x+y)+2
	TF2 * sinxpy = new TF2("sinxpy","sin(x+y)+2",xmin,xmax,ymin,ymax);
	benchmark -> SetTestFunction(sinxpy);

	// perform marginalization
	benchmark -> MCMCSetNIterationsRun(1000000);
	benchmark -> MarginalizeAll();

	// perform test
	benchmark -> PerformTest(
			benchmark -> GetBestFitParameters(),
			0,
			benchmark -> GetMarginalized("x","y"),
			true,
			"sinxpy.eps");




	// define and set test functions: gaus(x)*gaus(y) 
	TF2 * gaus2 = new TF2("gaus2","xygausn",xmin,xmax,ymin,ymax);
	gaus2->SetParameters(1,8,3,1,8,3);
	benchmark -> SetTestFunction(gaus2);

	// perform marginalization
	benchmark -> MCMCSetNIterationsRun(1000000);
	benchmark -> MarginalizeAll();

	// perform test
	benchmark -> PerformTest(
			benchmark -> GetBestFitParameters(),
			0,
			benchmark -> GetMarginalized("x","y"),
			true,
			"gaus2.eps");



	delete benchmark;

	return 0;

}

// ---------------------------------------------------------

