// RooFit headers needed before the other BAT headers (this is strange)
#include "RooRandom.h"

#include "BCRooInterface.h"

#include <BAT/BCLog.h>
#include <BAT/BCAux.h>
#include <BAT/BCH1D.h>
#include <BAT/BCH2D.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TRandom.h>

int main()
{
	// setting the RooFit random seed
	TRandom* myRandom = RooRandom::randomGenerator();
	myRandom->SetSeed(130);

	// set nice style for drawing than the ROOT default
	BCAux::SetStyle();

	// open log file with default level of logging
	BCLog::OpenLog("log.txt");
	BCLog::SetLogLevel(BCLog::detail);

	// create new RooFit-based model
	BCRooInterface * _myRooInterface = new BCRooInterface();
	_myRooInterface->Initialize("roodata.root");

	_myRooInterface -> MCMCSetNIterationsRun(1000000);

	// perform your analysis here
	_myRooInterface -> MarginalizeAll();
	_myRooInterface -> FindMode();
	_myRooInterface -> PrintAllMarginalized("plots.ps");
	_myRooInterface -> PrintResults("results.txt");

	std::cout << "\nThe results of the BAT calculations are:\n";
	std::cout << " Signal  " << std::endl;
	std::cout << " Mean  " << _myRooInterface -> GetMarginalized("S") -> GetMean() << std::endl;
	std::cout << " Median  " << _myRooInterface -> GetMarginalized("S") -> GetMedian() << std::endl;
	std::cout << " Mode  " << _myRooInterface -> GetMarginalized("S") -> GetMode() << std::endl;
	std::cout << " Mode  " << (_myRooInterface -> GetBestFitParameters()).at(0) << std::endl;
	std::cout << " Quantile 0.16  " << _myRooInterface -> GetMarginalized("S") -> GetQuantile(0.16) << std::endl;
	std::cout << " Quantile 0.84  " << _myRooInterface -> GetMarginalized("S") -> GetQuantile(0.84) << std::endl;
	std::cout << " Quantile 0.90  " << _myRooInterface -> GetMarginalized("S") -> GetQuantile(0.90) << std::endl;
	std::cout << " Quantile 0.95  " << _myRooInterface -> GetMarginalized("S") -> GetQuantile(0.95) << std::endl;
	std::cout << "\nCheck plots.ps and results.txt for more information on the results\n";

//	TFile * f = new TFile("plots.root", "RECREATE");
//	f -> cd();
//	_myRooInterface -> GetMarginalized("S") -> GetHistogram() -> Write();
//	_myRooInterface -> GetMarginalized("B") -> GetHistogram() -> Write();
//	_myRooInterface -> GetMarginalized("S", "B") -> GetHistogram() -> Write();
//	f -> Write();
//	f -> Close();

	// close log file
	BCLog::CloseLog();

	return 0;

}

