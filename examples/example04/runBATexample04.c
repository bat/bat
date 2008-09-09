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

	//	fModelEfficiency -> MCMCSetNIterationsRun(10000); 
	fModelEfficiency -> MCMCSetWriteChainToFile(true); 

	BCModelOutput * fModelEfficiencyOutput = new BCModelOutput(fModelEfficiency, "output_efficiency.root"); 

  // ---------------------------------------------------------
  // read data from file 
  // ---------------------------------------------------------

  BCDataSet * fDataSet = new BCDataSet(); 

  if (fDataSet -> ReadDataFromFileTxt("./data/data.txt", 2) != 0)
    return -1; 

  // assign data set to model 

  fModelEfficiency -> SetDataSet(fDataSet); 

  // ---------------------------------------------------------
  // optimize 
  // ---------------------------------------------------------

	fModelEfficiency -> SetOptimizationMethod(BCIntegrate::kOptMinuit); 

	fModelEfficiency -> FindMode(); 

  // ---------------------------------------------------------
  // marginalize 
  // ---------------------------------------------------------

	fModelEfficiency -> MarginalizeAll(); 

	//	fModelEfficiency -> GetMarginalized("epsilon") -> Print("modelefficiency_epsilon.ps",2);

	fModelEfficiency -> PrintAllMarginalized1D("marginalized"); 

  // ---------------------------------------------------------
  // evaluate p-value 
  // ---------------------------------------------------------

	double Nmax = fDataSet -> GetDataPoint(0) -> GetValue(1); 

	fModelEfficiency -> SetDataBoundaries(0, 0.0, Nmax); 
	fModelEfficiency -> FixDataAxis(1, true); 

	double loglikelihoodmax = fModelEfficiency -> LogLikelihood(fModelEfficiency -> GetBestFitParameters()); 

	fModelEfficiency -> CalculatePValue(fModelEfficiency -> GetBestFitParameters(), true) -> Print("pvalue.ps", 1, loglikelihoodmax); 
	
  // ---------------------------------------------------------
  // summarize
  // ---------------------------------------------------------

	fModelEfficiency -> PrintResults("summary.txt"); 

  // ---------------------------------------------------------
  // close log file 
  // ---------------------------------------------------------

	fModelEfficiencyOutput -> FillAnalysisTree(); 
	fModelEfficiencyOutput -> WriteMarginalizedDistributions(); 

	// write to file and close 

	fModelEfficiencyOutput -> Close(); 

  BCLog::CloseLog(); 

  return 0; 

}

// ---------------------------------------------------------
  
