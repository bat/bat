#include <BCModelEfficiency.h>
#include <BCModelOutput.h> 
#include <BCLog.h> 
#include <BCH1D.h> 

#include <TCanvas.h> 
#include <TGraphErrors.h>
#include <TH1D.h> 
#include <TH2D.h> 
#include <TF1.h> 

#include "style.c"

// ---------------------------------------------------------
  
int main()
{
  // ---------------------------------------------------------
  // open log file 
  // ---------------------------------------------------------

  BCLog::OpenLog("log.txt", BCLog::detail, BCLog::detail); 

  // ---------------------------------------------------------
  // set personal root style
  // ---------------------------------------------------------

  SetStyle();

  // ---------------------------------------------------------
  // model definition 
  // ---------------------------------------------------------

  BCModelEfficiency * fModelEfficiency = new BCModelEfficiency("ModelEfficiency"); 

	fModelEfficiency -> MCMCSetNChains(10); 
	//  fModelEfficiency -> MCMCSetFlagPCA(false); 
	//  fModelEfficiency -> MCMCSetNIterationsBurnIn(10000); 
	fModelEfficiency -> MCMCSetNIterationsRun(1000000); 
	//  fModelEfficiency -> MCMCSetTrialFunctionScale(0.01); 
	//  fModelEfficiency -> MCMCSetFlagInitialPosition(1); 
	//	fModelEfficiency -> MCMCSetWriteChainToFile(true); 

	BCModelOutput * fModelEfficiencyOutput = new BCModelOutput(fModelEfficiency, "output_efficiency.root"); 

	//	return 1; 

  // ---------------------------------------------------------
  // read data from file 
  // ---------------------------------------------------------

  BCDataSet * fDataSet = new BCDataSet(); 

  if (fDataSet -> ReadDataFromFileTxt("./data/data.txt", 2) != 0)
    return -1; 

  // assign data set to model 

  fModelEfficiency -> SetDataSet(fDataSet); 

	//	fModelEfficiency -> MCMCInitialize(); 

  // ---------------------------------------------------------
  // find maximum
  // ---------------------------------------------------------

	//  fModelEfficiency -> MCMCMetropolisHastings();  
	//	fModelEfficiency -> MCMCSimulatedAnnealing();  
  
	fModelEfficiency -> FindMode(); 

	fModelEfficiency -> MCMCMetropolis();  

  // ---------------------------------------------------------
  // marginalize 
  // ---------------------------------------------------------

	//  fModelEfficiency -> SetNbins(100);
	//  fModelEfficiency -> MarginalizeAll(); 

	fModelEfficiency -> GetMarginalized("epsilon") -> Print("modelefficiency_epsilon.ps",2);

  // ---------------------------------------------------------
  // summarize
  // ---------------------------------------------------------

	fModelEfficiency -> PrintSummary(); 

  // ---------------------------------------------------------
  // close log file 
  // ---------------------------------------------------------

	fModelEfficiencyOutput -> FillAnalysisTree(); 
	//	fModelEfficiencyOutput -> WriteMarginalizedDistributions(); 

	// write to file and close 

	fModelEfficiencyOutput -> Close(); 

  BCLog::CloseLog(); 

  return 0; 

}

// ---------------------------------------------------------
  
