#include <include/TestSuite.h> 
#include <include/PerfTest1DFunction.h>
#include <include/PerfTest2DFunction.h>

#include <TROOT.h>
#include <TMath.h>
#include <TF1.h>

#include <iostream>

using namespace std; 

int main()
{
	// create new test suite 
	TestSuite * ts = new TestSuite(); 

	//______________________________________________________________________________
	// add tests

	//______________________
	// 1D flat
	TF1* testfunc_1d_flat = new TF1("Flat", "1", -5., 5.0);
	PerfTest1DFunction*	perftest_1d_flat = new PerfTest1DFunction("1d_flat", testfunc_1d_flat); 
	perftest_1d_flat->GetSubtest("mode")->SetStatusOff(true);
	ts->AddTest(perftest_1d_flat); 

	//______________________
	// 1D slope
	TF1* testfunc_1d_slope = new TF1("Slope", "x", 0., 10.);
	PerfTest1DFunction*	perftest_1d_slope = new PerfTest1DFunction("1d_slope", testfunc_1d_slope); 
	ts->AddTest(perftest_1d_slope); 

	//______________________
	// 1D squared
	TF1* testfunc_1d_squared = new TF1("Squared", "400.-x*x", -20., 20.);
	PerfTest1DFunction*	perftest_1d_squared = new PerfTest1DFunction("1d_squared", testfunc_1d_squared); 
	ts->AddTest(perftest_1d_squared); 

	//______________________
	// 1D Gaussian
	TF1* testfunc_1d_gaus = new TF1("Gaus", "1.0/sqrt(2.0*TMath::Pi())/[1] * exp(-(x-[0])*(x-[0])/2/[1]/[1])", -25., 25.);
	testfunc_1d_gaus->FixParameter(0, 0.0);
	testfunc_1d_gaus->FixParameter(1, 5.0);
	PerfTest1DFunction*	perftest_1d_gaus = new PerfTest1DFunction("1d_gaus", testfunc_1d_gaus); 
	ts->AddTest(perftest_1d_gaus); 

	//______________________
	// 1D Poissons

	for (int i = 20; i <= 20; i++) {
		double xmax = 15; 
		if (i > 3)
			xmax = 10.0 * sqrt(double(i)); 
		TF1* testfunc_1d_poisson = new TF1("Poisson", "TMath::PoissonI([0], x)", 0., xmax);
		testfunc_1d_poisson->FixParameter(0, 10);
		PerfTest1DFunction* perftest = new PerfTest1DFunction(Form("1d_poisson_%i", i), testfunc_1d_poisson); 
		ts->AddTest(perftest); 
	}

	//______________________
	// 1D Binomials

	for (int N = 2; N < 3; N++) {
		for (int k = 0; k <= N; ++k) {
			TF1* testfunc = new TF1("Binomial", "([0]+1)*TMath::Binomial([0], [1]) * TMath::Power(x, [1]) * TMath::Power(1-x, [0]-[1])", 0., 1.);
			testfunc->FixParameter(0, N);
			testfunc->FixParameter(1, k);
			PerfTest1DFunction* perftest = new PerfTest1DFunction(Form("1d_binomial_%i_%i", k, N), testfunc); 
			ts->AddTest(perftest); 
		}
	}

	//______________________
	// 1D exponential

	TF1* testfunc_1d_exponential = new TF1("Exponential", "1/[0]*exp(-x/[0])", 0., 100.);
	testfunc_1d_exponential->FixParameter(0, 5);
	PerfTest1DFunction* perftest_1d_exponential = new PerfTest1DFunction("1d_exponential", testfunc_1d_exponential); 
	ts->AddTest(perftest_1d_exponential); 

	//______________________
	// 1D Cauchy

	TF1* testfunc_1d_cauchy = new TF1("Cauchy", "[1] / (3.14159 * ( (x-[2])**2 +[1]**2))", -25., 25.);
	testfunc_1d_cauchy->FixParameter(0, 0.);
	testfunc_1d_cauchy->FixParameter(1, 5.);
	PerfTest1DFunction* perftest_1d_cauchy = new PerfTest1DFunction("1d_cauchy", testfunc_1d_cauchy); 
	ts->AddTest(perftest_1d_cauchy); 

	//______________________
	// 1D Lognormal

	TF1* testfunc_1d_lognormal = new TF1("Lognormal", "1./sqrt(2*TMath::Pi()*[1])*1/x*exp(-(log(x)-[0])*(log(x)-[0])/2/[1]/[1])", 0., 10.);
	testfunc_1d_lognormal->FixParameter(0, 0.);
	testfunc_1d_lognormal->FixParameter(1, 1.);
	PerfTest1DFunction* perftest_1d_lognormal = new PerfTest1DFunction("1d_lognormal", testfunc_1d_lognormal); 
	ts->AddTest(perftest_1d_lognormal); 

	//______________________
	// 1D x^4 sin^2(x)

	TF1* testfunc_1d_sin2 = new TF1("Sin2", "x*x*x*x*sin(x)*sin(x)", 2., 25.);
	PerfTest1DFunction* perftest_1d_sin2 = new PerfTest1DFunction("1d_sin2", testfunc_1d_sin2); 
	ts->AddTest(perftest_1d_sin2); 

	//______________________
	// 1D 2 Gaussians
	TF1* testfunc_1d_2gaus = new TF1("2gaus1d", "gaus + gaus(3)", -25., 50.);
	testfunc_1d_2gaus->FixParameter(0,  1.0);
	testfunc_1d_2gaus->FixParameter(1, -10.0);
	testfunc_1d_2gaus->FixParameter(2,  2.0);
	testfunc_1d_2gaus->FixParameter(3,  2.0);
	testfunc_1d_2gaus->FixParameter(4, 30.0);
	testfunc_1d_2gaus->FixParameter(5,  1.0);
	PerfTest1DFunction*	perftest_1d_2gaus = new PerfTest1DFunction("1d_2gaus", testfunc_1d_2gaus); 
	perftest_1d_2gaus->SetNbins("x", 200);
	ts->AddTest(perftest_1d_2gaus); 

	//______________________
	// 2D flat
 	TF2* testfunc_2d_flat = new TF2("Flat", "1", -5., 5., -5., 5.);
	PerfTest2DFunction*	perftest_2d_flat = new PerfTest2DFunction("2d_flat", testfunc_2d_flat); 
	ts->AddTest(perftest_2d_flat); 

	//______________________
	// 2D Gaussian
	TF2* testfunc_2d_gaus = new TF2("Gaus", "xygausn", -3., 3., -5., 7.);
	testfunc_2d_gaus->SetParameters(1,0,1,1,2);
	PerfTest2DFunction*	perftest_2d_gaus = new PerfTest2DFunction("2d_gaus", testfunc_2d_gaus); 
	ts->AddTest(perftest_2d_gaus); 

	//______________________
	// 1D 2 Gaussians
	TF2* testfunc_2d_2gaus = new TF2("2twoGaus2d", 
																	 "[0] * ( [1]*exp(-0.5*((x-[2])/[3])**2)*exp(-0.5*((y-[4])/[5])**2) + [6]*exp(-0.5*((x-[7])/[8])**2)*exp(-0.5*((y-[9])/[10])**2))",
																	 -20., 20., -20., 20);
	testfunc_2d_2gaus->SetParameters(1.,   10., 0., 1.0,  5., 1.0,    10., 5., 1.0,  10., 1.0);
	PerfTest2DFunction*	perftest_2d_2gaus = new PerfTest2DFunction("2d_2gaus", testfunc_2d_2gaus); 
	ts->AddTest(perftest_2d_2gaus); 

	//______________________________________________________________________________
	// perform all tests 
	ts -> RunTests(); 

	// print results to screen
	ts -> PrintResultsScreen(); 

	// print results to html
	ts -> PrintResultsHTML(); 

	// delete test suite 
	delete ts;  

	return 0; 
}
