#include <include/ReleaseTestSuite.h> 
#include <include/PerfTest.h>

#include <TROOT.h>
#include <TMath.h>
#include <TF1.h>

#include <iostream>

using namespace std; 

int main()
{
	// create new test suite 
	ReleaseTestSuite* rts = new ReleaseTestSuite(); 

	// prepare test suite
	rts->PrepareTests(); 

	// set precision
	//	rts->SetPrecision(PerfTest::kDetail);
	rts->SetPrecision(PerfTest::kCoarse);

	// run all tests
	rts->RunTests();

	// print results to screen
	rts->PrintResultsScreen(); 

	// print results to html
	rts->PrintResultsHTML(); 

	// delete test suite 
	delete rts;  

	return 0; 
}
