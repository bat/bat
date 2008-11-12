#include "BCLog.h"
#include "BCBenchmarkMCMC.h"
#include "style.c"

#include <TF1.h>

// ---------------------------------------------------------

int main()
{
	SetStyle();

	BCLog::OpenLog("log0.txt", BCLog::detail, BCLog::detail);

	BCBenchmarkMCMC * benchmark = new BCBenchmarkMCMC("benchmark");

	// add parameters
	double xmin =  0.0;
	double xmax = 50.0;
	benchmark -> AddParameter("x",  xmin, xmax);
	benchmark -> SetNbins(500, 0);

	// define and set test functions
//	TF1 * testfunction1 = new TF1("testfunction1", "x*x", xmin, xmax);
	TF1 * testfunction2 = new TF1("testfunction2", "x*x*sin(x)*sin(x)", xmin, xmax);
	benchmark -> SetTestFunction(testfunction2);

	// perform marginalization
	benchmark -> MCMCSetNIterationsRun(1000000);
	benchmark -> MarginalizeAll();

	// perform test
	benchmark -> PerformTest(
			benchmark -> GetBestFitParameters(),
			0,
			benchmark -> GetMarginalized("x"),
			true,
			"test.ps");

	delete benchmark;

	return 0;
}

// ---------------------------------------------------------

