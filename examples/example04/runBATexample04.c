#include <BCModelEfficiency.h>
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

  BCLog::OpenLog(); 

  // ---------------------------------------------------------
  // set personal root style
  // ---------------------------------------------------------

  SetStyle();

  // ---------------------------------------------------------
  // model definition 
  // ---------------------------------------------------------

  BCModelEfficiency * fModelEfficiency = new BCModelEfficiency("ModelEfficiency"); 

  // ---------------------------------------------------------
  // read data from file 
  // ---------------------------------------------------------

  BCDataSet * fDataSet = new BCDataSet(); 

  if (fDataSet -> ReadDataFromFileTxt("./data/data.txt", 2) != 0)
    return -1; 

  // assign data set to model 

  fModelEfficiency -> SetDataSet(fDataSet); 

  // ---------------------------------------------------------
  // find maximum
  // ---------------------------------------------------------

  fModelEfficiency -> FindMode(); 

  // ---------------------------------------------------------
  // marginalize 
  // ---------------------------------------------------------

  fModelEfficiency -> SetNbins(100);
  fModelEfficiency -> MarginalizeAll(); 

  fModelEfficiency -> GetMarginalized("efficiency") -> Print("modelefficiency_epsilon.ps");
  fModelEfficiency -> GetMarginalized("efficiency") -> Print("modelefficiency_epsilon_1.ps",2);

  // ---------------------------------------------------------
  // summarize
  // ---------------------------------------------------------

  fModelEfficiency -> PrintSummary(); 

  // ---------------------------------------------------------
  // close log file 
  // ---------------------------------------------------------

  BCLog::CloseLog(); 

  return 0; 

}

// ---------------------------------------------------------
  
