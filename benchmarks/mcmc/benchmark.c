#include <BCLog.h>

#include "style.c"
#include <BCBenchmarkMCMC.h> 
#include "TF1.h" 

// ---------------------------------------------------------
  
int main()
{

	// ---------------------------------------------------------
	// set style  
	// ----------------------------------------------------------

	// calls a function which defines a nicer style than the ROOT
	// default.
	SetStyle(); 

	// ---------------------------------------------------------
	// open log file 
	// ---------------------------------------------------------

	// opens the log file. 
	BCLog::OpenLog("log.txt", BCLog::detail, BCLog::detail); 

	// ---------------------------------------------------------
	// create model 
	// ---------------------------------------------------------

	BCBenchmarkMCMC * benchmark = new BCBenchmarkMCMC("benchmark"); 

	// add parameters 
	double xmin =  0.0; 
	double xmax = 50.0; 
	benchmark -> AddParameter("x",  xmin, xmax);
	benchmark -> SetNbins(500, 0); 

	// define test functions 
	TF1 * testfunction1 = new TF1("testfunction1", "x*x", xmin, xmax); 
	TF1 * testfunction2 = new TF1("testfunction2", "x*x*sin(x)*sin(x)", xmin, xmax); 

	// set test function 
	benchmark -> SetTestFunction(testfunction2); 

	// perform marginalization 
	benchmark -> MarginalizeAll(); 

	BCH1D * hist = benchmark -> GetMarginalized("x"); 
	
	// perform test 
	double* chi2; 
	benchmark -> PerformTest(benchmark -> GetBestFitParameters(), 
													 0, 
													 hist,
													 chi2, 
													 true, 
													 "test.ps"); 

	delete benchmark; 

}

// ---------------------------------------------------------
  
