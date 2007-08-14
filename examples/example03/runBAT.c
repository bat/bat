#include <BCModelBackground.h>
#include <BCModelSignal.h>
#include <BCParameter.h>
#include <BCModelManager.h>
#include <BCDataPoint.h> 
#include <BCH1D.h> 

#include <TCanvas.h> 
#include <TH1D.h> 
#include <TF1.h> 
#include <BCLog.h> 

#include  <iostream>

#include "style.c"

// ---------------------------------------------------------
  
int main()
{

  // open log file 
  
  BCLog::OpenLog(); 

  // set root style 

  SetStyle(); 

  // ---------------------------------------------------------
  // define models 
  // ---------------------------------------------------------

  BCModelBackground* fModelBackground = new BCModelBackground("background model"); 
  BCModelSignal* fModelSignal = new BCModelSignal("signal model"); 

  // ---------------------------------------------------------
  // define manager 
  // ---------------------------------------------------------

  BCModelManager* fModelManager = new BCModelManager(); 

  // ---------------------------------------------------------
  // add models to model manager with a priori probabilities 
  // ---------------------------------------------------------

  fModelManager -> AddModel(fModelBackground, 0.5); 
  fModelManager -> AddModel(fModelSignal, 0.5); 

  // ---------------------------------------------------------
  // read data from file 
  // ---------------------------------------------------------

  BCDataSet* fDataSet = new BCDataSet(); 

  // read data from ROOT file 
  
  int errorcode = fDataSet -> ReadDataFromFileTxt("./data/data.txt", 2); 

  if (errorcode != 0)
    return -1;
  
  // connect mananger (and models) with data set
  
  fModelManager -> SetDataSet(fDataSet);

  // ---------------------------------------------------------
  // initialize
  // ---------------------------------------------------------

  fModelManager -> Initialize();

  // ---------------------------------------------------------
  // set boundaries on possible data values  
  // ---------------------------------------------------------

  BCDataPoint* fDataPointLowerBoundaries = new BCDataPoint(2); 
  fDataPointLowerBoundaries -> SetValue(0, fDataSet -> GetDataPoint(0) -> GetValue(0)); 
  fDataPointLowerBoundaries -> SetValue(1, 0.0); 

  BCDataPoint* fDataPointUpperBoundaries = new BCDataPoint(2); 
  fDataPointUpperBoundaries -> SetValue(0, fDataSet -> GetDataPoint(0) -> GetValue(1)); 
  fDataPointUpperBoundaries -> SetValue(1, 10.0); 

  int nmimimum = fDataSet -> GetNDataPoints() - 1; 
  int nmaximum = fDataSet -> GetNDataPoints() - 1; 

  // set boundaries in models 

  fModelBackground -> SetDataPointLowerBoundaries(fDataPointLowerBoundaries); 
  fModelBackground -> SetDataPointUpperBoundaries(fDataPointUpperBoundaries); 
  fModelBackground -> SetNDataPointsMinimum(nmimimum); 
  fModelBackground -> SetNDataPointsMaximum(nmaximum); 

  fModelSignal -> SetDataPointLowerBoundaries(fDataPointLowerBoundaries); 
  fModelSignal -> SetDataPointUpperBoundaries(fDataPointUpperBoundaries); 
  fModelSignal -> SetNDataPointsMinimum(nmimimum); 
  fModelSignal -> SetNDataPointsMaximum(nmaximum); 

  // ---------------------------------------------------------
  // marginalize 
  // ---------------------------------------------------------

  fModelBackground -> MarginalizeAll();
  fModelBackground -> GetMarginalized("background") -> Print("modelbackground_background.ps");

  fModelSignal -> MarginalizeAll();
  fModelSignal -> GetMarginalized("background") -> Print("modelsignal_background.ps");
  fModelSignal -> GetMarginalized("signal") -> Print("modelsignal_signal.ps");
  fModelSignal -> GetMarginalized("background", "signal") -> Print("modelsignal_background_signal.ps", 2);

  // ---------------------------------------------------------
  // summarize
  // ---------------------------------------------------------
  
  fModelManager -> PrintSummary(); 

  fModelManager -> PrintSummaryToFile("output.log"); 

  // ---------------------------------------------------------
  // Print data with best fit result 
  // ---------------------------------------------------------


  // close log file 

  BCLog::CloseLog(); 

  return 0; 

}

// ---------------------------------------------------------
  
