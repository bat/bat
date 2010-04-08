// RooFit headers needed before the other BAT headers (this is strange)
#include "RooRandom.h"
#include "RooWorkspace.h"

#include "BCRooInterface.h"

#include <BAT/BCLog.h>
#include <BAT/BCAux.h>
#include <BAT/BCH1D.h>
#include <BAT/BCH2D.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TRandom.h>

int main(int argc, char **argv)
{
  std::cout << "Running BAT via the BCRooInterface...\n";

  // retrieve the command line arguments
  char* rootFile = "bat_workspace.root";
  char* wsName = "batWS";
  char* outputFile = "bat_plots.ps";
  int nMCMC = 1000000;

  // display a quick help if no arguments are specified
  if (argc==0) {
    std::cout << "Run with: ./runRooInterface bat_workspace.root batWS bat_plots.ps\n";
    std::cout << "the four arguments are:\n";
    std::cout << " - the name of the ROOT file that contains a workspace. This workspace should contain:\n";
    std::cout << "  - data: a RooAbsData object holding the data measurements\n";
    std::cout << "  - model: a RooAbsPdf object holding the probability density function\n";
    std::cout << "  - priorPOI and priorNuisance: two RooAbsPdf objects representing respectively the prior probability on the parameter of interest and the nuisance parameters.\n";
    std::cout << "  - POI: a RooArgList holding one single object being the parameter of interest\n";
    std::cout << "  - parameters: a RooArgList holding all nuisance parameters\n";
    std::cout << " - the name workspace to retrieve from the ROOT file\n";
    std::cout << " - the name of the postscript file that will contain the resulting posterior probability plots\n";
    std::cout << " - the number of MCMC iterations (default = 1000000)\n";
    std::cout << "For more information see the README file of the BCRooInterface\n";
  }

  if (argc>=1) rootFile = argv[1];
  if (argc>=2) wsName = argv[2];
  //if (argc>=3) outputFile = argv[3];
  //if (argc>=4) nMCMC = int(argv[4]);

  std::cout << "The inputs will be retrieved from " << rootFile << " (workspace " << wsName << ") and the posterior probability plot will be stored in " << outputFile << std::endl;
  std::cout << nMCMC << " MCMC iteration will be performed\n";

  // set nice style for drawing than the ROOT default
  BCAux::SetStyle();

  // open log file with default level of logging
  BCLog::OpenLog("bat_log.txt");
  BCLog::SetLogLevel(BCLog::detail);

  // create new RooFit-based model
  BCRooInterface * _myRooInterface = new BCRooInterface();
  _myRooInterface->Initialize(rootFile,wsName,"data","model","priorPOI","priorNuisance","parameters","POI");

  _myRooInterface -> MCMCSetNIterationsRun(nMCMC);

  // perform your analysis here
  _myRooInterface -> MarginalizeAll();
  _myRooInterface -> FindMode();
  _myRooInterface -> PrintAllMarginalized(outputFile);
  _myRooInterface -> PrintResults("bat_results.txt");

  TFile* file = new TFile(rootFile);
  RooWorkspace* bat_ws = (RooWorkspace*) file->Get(wsName);
  TString parameterName = bat_ws->set("POI")->first()->GetName();

  std::cout << "\nThe results of the BAT calculations are:\n";
  std::cout << " Signal  " << std::endl;
  std::cout << " Mean  " << _myRooInterface -> GetMarginalized(parameterName) -> GetMean() << std::endl;
  std::cout << " Median  " << _myRooInterface -> GetMarginalized(parameterName) -> GetMedian() << std::endl;
  std::cout << " Mode  " << _myRooInterface -> GetMarginalized(parameterName) -> GetMode() << std::endl;
  std::cout << " Mode  " << (_myRooInterface -> GetBestFitParameters()).at(0) << std::endl;
  std::cout << " Quantile 0.16  " << _myRooInterface -> GetMarginalized(parameterName) -> GetQuantile(0.16) << std::endl;
  std::cout << " Quantile 0.84  " << _myRooInterface -> GetMarginalized(parameterName) -> GetQuantile(0.84) << std::endl;
  std::cout << " Quantile 0.90  " << _myRooInterface -> GetMarginalized(parameterName) -> GetQuantile(0.90) << std::endl;
  std::cout << " Quantile 0.95  " << _myRooInterface -> GetMarginalized(parameterName) -> GetQuantile(0.95) << std::endl;
  std::cout << "\nCheck " << outputFile << " and results.txt for more information on the results\n";

  // close log file
  BCLog::CloseLog();

  return 0;

}

