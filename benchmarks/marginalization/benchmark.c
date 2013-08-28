#include <BCModelOutput.h>
#include <BCModelGauss.h>
#include <BCLog.h>

#include <TStopwatch.h>

// ---------------------------------------------------------

int main()
{

	// ---------------------------------------------------------
	// create stopwatch
	// ---------------------------------------------------------

	TStopwatch * fStopwatch = new TStopwatch();

	// ---------------------------------------------------------
	// open log file
	// ---------------------------------------------------------

	BCLog::OpenLog();

	// ---------------------------------------------------------
	// model definition
	// ---------------------------------------------------------

	BCModelGauss* fModelGauss = new BCModelGauss("ModelGauss");

	BCModelOutput * fModelOutputGauss = new BCModelOutput(fModelGauss, "output.root");
	fModelOutputGauss -> WriteMarkovChain(true);

	// ---------------------------------------------------------
	// do marginalization with different methods
	// ---------------------------------------------------------

	// Monte Carlo integration

	// start stopwatch

	std::cout << std::endl;
	std::cout << " Marginalize with Markov chain " << endl;

	fModelGauss -> SetMarginalizationMethod(BCIntegrate::kMMetropolis);
	fModelGauss -> SetNbins(100);

	fStopwatch -> Start();

	fModelGauss -> MarginalizeAll();

	fStopwatch -> Stop();

	double realtime_markov = fStopwatch -> RealTime();
	double cputime_markov  = fStopwatch -> CpuTime();

	fModelGauss -> GetMarginalized("x1") -> Print("x1.pdf");

	// ---------------------------------------------------------
	// write to output file
	// ---------------------------------------------------------

	// fill the ROOT file with the actual output of the model.

	fModelOutputGauss -> FillAnalysisTree();

	fModelOutputGauss -> WriteMarginalizedDistributions();

	// write to file and close

	fModelOutputGauss -> Close();

	// other integration methods ...

	// ..

	// ---------------------------------------------------------
	// print results of benchmark test
	// ---------------------------------------------------------

	std::cout << std::endl;
	std::cout << " Summary of integration methods " << std::endl;
	std::cout << " ============================== " << std::endl;
	std::cout << std::endl;

	std::cout << " Markov chains " << std::endl;
	std::cout << " Real time [s] : " << realtime_markov << std::endl;
	std::cout << " CPU time  [s] : " << cputime_markov << std::endl;
	std::cout << std::endl;

	// ---------------------------------------------------------
	// close log file
	// ---------------------------------------------------------

	// closes the log file

	BCLog::CloseLog();

	return 0;

}

// ---------------------------------------------------------

