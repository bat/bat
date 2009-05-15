// test MCMC sampling for 1D functions
// Authors: Kevin Kroeninger, Jing Liu
//
#include "BAT/BCAux.h"
#include "BAT/BCLog.h"

#include <TF1.h>

#include "BCBenchmarkMCMC.h"

// ---------------------------------------------------------

int main()
{
	BCAux::SetStyle();

	BCLog::OpenLog("SamplingTest1D.log", BCLog::detail, BCLog::detail);

	BCBenchmarkMCMC * benchmark = new BCBenchmarkMCMC("benchmark");

	// add parameters
	double xmin = -4.0;
	double xmax = 16.0;
	benchmark -> AddParameter("x", xmin, xmax);
	//benchmark -> SetNbins(500, 0);




	// define and set test functions: pol0
	TF1 * testpol0 = new TF1("Constant", "pol0(0)", xmin, xmax);
	testpol0->SetParameter(0,6);
	benchmark -> SetTestFunction(testpol0);

	// perform marginalization
	benchmark -> MCMCSetNIterationsRun(1000000);
	benchmark -> MarginalizeAll();

	// perform test
	benchmark -> PerformTest(
			benchmark -> GetBestFitParameters(),
			0,
			benchmark -> GetMarginalized("x"),
			true,
			"pol0.eps");




	// define and set test functions: pol1
	TF1 * testpol1 = new TF1("1-order Polynomial", "pol1(0)", xmin, xmax);
	testpol1->SetParameters(2,3);
	benchmark -> SetTestFunction(testpol1);

	// perform marginalization
	benchmark -> MCMCSetNIterationsRun(1000000);
	benchmark -> MarginalizeAll();

	// perform test
	benchmark -> PerformTest(
			benchmark -> GetBestFitParameters(),
			0,
			benchmark -> GetMarginalized("x"),
			true,
			"pol1.eps");





	// define and set test functions: pol2
	TF1 * testpol2 = new TF1("2-order Polynomial", "pol2(0)", xmin, xmax);
	testpol2->SetParameters(2,3,5);
	benchmark -> SetTestFunction(testpol2);

	// perform marginalization
	benchmark -> MCMCSetNIterationsRun(1000000);
	benchmark -> MarginalizeAll();

	// perform test
	benchmark -> PerformTest(
			benchmark -> GetBestFitParameters(),
			0,
			benchmark -> GetMarginalized("x"),
			true,
			"pol2.eps");





	// define and set test functions: exponential decay
	TF1 * testexp = new TF1("Exponential decay", "[0]*exp(-[1]*x)", xmin, xmax);
	testexp->SetParameters(1,0.2);
	benchmark -> SetTestFunction(testexp);

	// perform marginalization
	benchmark -> MCMCSetNIterationsRun(1000000);
	benchmark -> MarginalizeAll();

	// perform test
	benchmark -> PerformTest(
			benchmark -> GetBestFitParameters(),
			0,
			benchmark -> GetMarginalized("x"),
			true,
			"expdecay.eps");





	// define and set test functions: gaus
	TF1 * testgaus = new TF1("Gaussian", "gaus(0)", xmin, xmax);
	testgaus->SetParameters(2,3,4);
	benchmark -> SetTestFunction(testgaus);

	// perform marginalization
	benchmark -> MCMCSetNIterationsRun(1000000);
	benchmark -> MarginalizeAll();

	// perform test
	benchmark -> PerformTest(
			benchmark -> GetBestFitParameters(),
			0,
			benchmark -> GetMarginalized("x"),
			true,
			"gaus.eps");

	delete benchmark;

	return 0;
}

// ---------------------------------------------------------

