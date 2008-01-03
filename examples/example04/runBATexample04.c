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

  fModelEfficiency -> MCMCSetNChains(1); 
  fModelEfficiency -> MCMCSetFlagPCA(false); 
  fModelEfficiency -> MCMCSetNIterationsBurnIn(10000); 
  fModelEfficiency -> MCMCSetNIterationsMax(100000); 
	//  fModelEfficiency -> MCMCSetTrialFunctionScale(0.01); 
  fModelEfficiency -> MCMCSetFlagInitialPosition(1); 
	fModelEfficiency -> MCMCSetWriteChainToFile(true); 

	//	fModelEfficiency -> MCMCInitialize(); 

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

	

	fModelEfficiency -> MCMCMetropolis();  
	//  fModelEfficiency -> MCMCMetropolisHastings();  
	//	fModelEfficiency -> MCMCSimulatedAnnealing();  

  TCanvas * can1 = new TCanvas(); 
  can1 -> Divide(2, 1); 
  can1 -> cd(1); 
  fModelEfficiency -> hist_efficiency -> Draw(); 
  can1 -> cd(2); 
  fModelEfficiency -> hist_bestfit -> Draw(); 
  can1 -> Print("canvas_results.ps"); 
  
  //  fModelEfficiency -> FindMode(); 

  // ---------------------------------------------------------
  // marginalize 
  // ---------------------------------------------------------

	//  fModelEfficiency -> SetNbins(100);
	//  fModelEfficiency -> MarginalizeAll(); 

	fModelEfficiency -> GetMarginalized("mu1") -> Print("modelefficiency_epsilon.ps");
	fModelEfficiency -> GetMarginalized("mu1", "mu2") -> Print("modelefficiency_mu1_mu2.ps");
	//  fModelEfficiency -> GetMarginalized("efficiency") -> Print("modelefficiency_epsilon_1.ps",2);

  // ---------------------------------------------------------
  // summarize
  // ---------------------------------------------------------

	//  fModelEfficiency -> PrintSummary(); 

  // ---------------------------------------------------------
  // close log file 
  // ---------------------------------------------------------

	//	fModelEfficiencyOutput -> FillAnalysisTree(); 
	//	fModelEfficiencyOutput -> WriteMarginalizedDistributions(); 

	// write to file and close 

	fModelEfficiencyOutput -> Close(); 

  BCLog::CloseLog(); 

  return 0; 

}

// ---------------------------------------------------------
  
