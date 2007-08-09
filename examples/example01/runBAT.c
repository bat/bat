#include <BCModelPol1.h>
#include <BCLog.h> 

#include <TCanvas.h> 
#include <TGraphErrors.h>
#include <TH1D.h> 
#include <TH2D.h> 
#include <TF1.h> 

// ---------------------------------------------------------
  
int main()
{
  // ---------------------------------------------------------
  // open log file 
  // ---------------------------------------------------------

  BCLog::OpenLog(); 

  // ---------------------------------------------------------
  // model definition 
  // ---------------------------------------------------------

  BCModelPol1* fModelPol1 = new BCModelPol1("ModelPol1"); 

  // ---------------------------------------------------------
  // read data from file 
  // ---------------------------------------------------------

  BCDataSet* fDataSet = new BCDataSet(); 

  if (fDataSet -> ReadDataFromFileTxt("./data/data.txt", 3) != 0)
    return -1; 

  // assign data set to model 

  fModelPol1 -> SetDataSet(fDataSet); 

  // ---------------------------------------------------------
  // normalize  
  // ---------------------------------------------------------

//  fModelPol1 -> FindMode(); 

  // ---------------------------------------------------------
  // marginalize 
  // ---------------------------------------------------------

//  fModelPol1 -> SetNbins(200);
  fModelPol1 -> MarginalizeAll();

  fModelPol1 -> GetMarginalized("constant") -> Print("modelpol1_constant.ps", 3);
  fModelPol1 -> GetMarginalized("slope") -> Print("modelpol1_slope.ps", 3);
  fModelPol1 -> GetMarginalized("constant", "slope") -> Print("modelpol1_constant_slope.ps", 2);

  // ---------------------------------------------------------
  // summarize
  // ---------------------------------------------------------

  fModelPol1 -> PrintSummary(); 

  // ---------------------------------------------------------
  // Print data with best fit result 
  // ---------------------------------------------------------

  // canvas 

  TCanvas* canvas_bestfit = new TCanvas(); 
  canvas_bestfit -> cd(); 

  // axes histogram 

  TH2D* hist_axes = new TH2D("hist_axes", "Data;x;y", 1, 0.0, 100.0, 1, 0.0, 6.0); 
  hist_axes -> SetStats(false); 
  hist_axes -> Draw(); 

  // graph with data 

  TGraphErrors* graph = new TGraphErrors(); 
  graph -> SetMarkerStyle(20); 

  // set data points 

  for (int i = 0; i < fModelPol1 -> GetNDataPoints(); i++)
    {
      graph -> SetPoint(i, 
			fModelPol1 -> GetDataPoint(i) -> GetValue(0), 
			fModelPol1 -> GetDataPoint(i) -> GetValue(1)); 
      graph -> SetPointError(i, 
			     0.0, 
			     fModelPol1 -> GetDataPoint(i) -> GetValue(2)); 
    }

  if (fModelPol1 -> GetBestFitParameters().size() > 0)
    {
      // best fit function 
      
      TF1* func_pol1 = new TF1("func_pol1", "[0] + [1] * x", 0, 100); 
      
      func_pol1 -> SetParameter(0, fModelPol1 -> GetBestFitParameter(0)); 
      func_pol1 -> SetParameter(1, fModelPol1 -> GetBestFitParameter(1)); 

      func_pol1 -> Draw("SAME"); 
    }

  // draw graph 

  graph -> Draw("SAMEP"); 

  // print to file 

  canvas_bestfit -> Print("data.ps"); 

  // ---------------------------------------------------------
  // close log file 
  // ---------------------------------------------------------

  BCLog::CloseLog(); 

  return 0; 

}

// ---------------------------------------------------------
  
