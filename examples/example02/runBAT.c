#include <BCModelPol0.h>
#include <BCModelPol1.h>
#include <BCModelPol2.h>
#include <BCModelManager.h>
#include <BCDataPoint.h> 
#include <BCLog.h> 

#include <TCanvas.h> 
#include <TGraphErrors.h>
#include <TH1D.h> 
#include <TH2D.h> 
#include <TF1.h> 
#include <TLegend.h> 
#include <TLine.h> 
#include <TText.h> 

#include <iostream>

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

  BCModelPol0* fModelPol0 = new BCModelPol0("ModelPol0"); 
  BCModelPol1* fModelPol1 = new BCModelPol1("ModelPol1"); 
  BCModelPol2* fModelPol2 = new BCModelPol2("ModelPol2"); 

  // ---------------------------------------------------------
  // define model manager 
  // ---------------------------------------------------------

  BCModelManager* fModelManager = new BCModelManager(); 

  // ---------------------------------------------------------
  // add models to model manager with a priori probabilities 
  // ---------------------------------------------------------

  fModelManager -> AddModel(fModelPol0, 1.0/3.0); 
  fModelManager -> AddModel(fModelPol1, 1.0/3.0); 
  fModelManager -> AddModel(fModelPol2, 1.0/3.0); 

  // ---------------------------------------------------------
  // read data from file 
  // ---------------------------------------------------------

  BCDataSet* fDataSet = new BCDataSet(); 

  int errorcode = fDataSet -> ReadDataFromFileTxt("./data/data_ModelPol2.txt", 3); 

  if (errorcode != 0) 
    return -1; 

  // connect mananger (and models) with data set 
  
  fModelManager -> SetDataSet(fDataSet); 

  // ---------------------------------------------------------
  // initialize 
  // ---------------------------------------------------------

  fModelManager -> Initialize(); 

  // ---------------------------------------------------------
  // marginalize 
  // ---------------------------------------------------------

  if (fModelPol0 -> GetModelAPosterioriProbability() > fModelPol1 -> GetModelAPosterioriProbability() && 
      fModelPol0 -> GetModelAPosterioriProbability() > fModelPol2 -> GetModelAPosterioriProbability()) 
    fModelPol0 -> MarginalizeProbability("constant") -> Print("modelpol0_constant.ps", 1); 

  if (fModelPol1 -> GetModelAPosterioriProbability() > fModelPol0 -> GetModelAPosterioriProbability() && 
      fModelPol1 -> GetModelAPosterioriProbability() > fModelPol2 -> GetModelAPosterioriProbability()) 
    {
      fModelPol1 -> MarginalizeProbability("constant")          -> Print("modelpol1_constant.ps", 1); 
      fModelPol1 -> MarginalizeProbability("slope")             -> Print("modelpol1_slope.ps", 1); 
      fModelPol1 -> MarginalizeProbability("constant", "slope") -> Print("modelpol1_constant_slope.ps", 2); 
    }

  if (fModelPol2 -> GetModelAPosterioriProbability() > fModelPol0 -> GetModelAPosterioriProbability() && 
      fModelPol2 -> GetModelAPosterioriProbability() > fModelPol1 -> GetModelAPosterioriProbability()) 
    {
      fModelPol2 -> MarginalizeProbability("constant")          -> Print("modelpol2_constant.ps", 1); 
      fModelPol2 -> MarginalizeProbability("slope")             -> Print("modelpol2_slope.ps", 1); 
      fModelPol2 -> MarginalizeProbability("quad")              -> Print("modelpol2_quad.ps", 1); 
      fModelPol2 -> MarginalizeProbability("constant", "slope") -> Print("modelpol2_constant_slope.ps", 2); 
      fModelPol2 -> MarginalizeProbability("constant", "quad")  -> Print("modelpol2_constant_quad.ps", 2); 
      fModelPol2 -> MarginalizeProbability("slope", "quad")     -> Print("modelpol2_slope_quad.ps", 2); 
    }

  // ---------------------------------------------------------
  // Do goodness-of-fit test  
  // ---------------------------------------------------------

  // set boundaries on possible data values  

  BCDataPoint* fDataPointLowerBoundaries = new BCDataPoint(3); 
  fDataPointLowerBoundaries -> SetValue(0, 0.0); 
  fDataPointLowerBoundaries -> SetValue(1, 0.0); 
  fDataPointLowerBoundaries -> SetValue(2, 0.2); 

  BCDataPoint* fDataPointUpperBoundaries = new BCDataPoint(3); 
  fDataPointUpperBoundaries -> SetValue(0, 100.0); 
  fDataPointUpperBoundaries -> SetValue(1,  10.0); 
  fDataPointUpperBoundaries -> SetValue(2,   0.2); 

  // set boundaries in models 
  
  fModelPol0 -> SetDataPointLowerBoundaries(fDataPointLowerBoundaries); 
  fModelPol0 -> SetDataPointUpperBoundaries(fDataPointUpperBoundaries); 
  fModelPol0 -> SetNDataPointsMinimum(10); 
  fModelPol0 -> SetNDataPointsMaximum(10); 

  fModelPol1 -> SetDataPointLowerBoundaries(fDataPointLowerBoundaries); 
  fModelPol1 -> SetDataPointUpperBoundaries(fDataPointUpperBoundaries); 
  fModelPol1 -> SetNDataPointsMinimum(10); 
  fModelPol1 -> SetNDataPointsMaximum(10); 

  fModelPol2 -> SetDataPointLowerBoundaries(fDataPointLowerBoundaries); 
  fModelPol2 -> SetDataPointUpperBoundaries(fDataPointUpperBoundaries); 
  fModelPol2 -> SetNDataPointsMinimum(10); 
  fModelPol2 -> SetNDataPointsMaximum(10); 

 // set up grid

  std::vector <bool> grid;
  grid.push_back(true);
  grid.push_back(false);
  grid.push_back(false);

  std::vector <double> limits;
  limits.push_back( 5.0);
  limits.push_back(10.0);
  limits.push_back(10.0);

  if (fModelPol0 -> GetModelAPosterioriProbability() > fModelPol1 -> GetModelAPosterioriProbability() && 
      fModelPol0 -> GetModelAPosterioriProbability() > fModelPol2 -> GetModelAPosterioriProbability()) 
    fModelPol0 -> DoGoodnessOfFitTest(100, fModelPol0 -> GetBestFitParameters(), grid, limits) -> 
      Print("modelpol0_gof.ps", 2, TMath::Log10(fModelPol0 -> Likelihood(fModelPol0 -> GetBestFitParameters())));

  if (fModelPol1 -> GetModelAPosterioriProbability() > fModelPol0 -> GetModelAPosterioriProbability() && 
      fModelPol1 -> GetModelAPosterioriProbability() > fModelPol2 -> GetModelAPosterioriProbability()) 
    fModelPol1 -> DoGoodnessOfFitTest(100, fModelPol1 -> GetBestFitParameters(), grid, limits) -> 
      Print("modelpol1_gof.ps", 2, TMath::Log10(fModelPol1 -> Likelihood(fModelPol1 -> GetBestFitParameters())));

  if (fModelPol2 -> GetModelAPosterioriProbability() > fModelPol0 -> GetModelAPosterioriProbability() && 
      fModelPol2 -> GetModelAPosterioriProbability() > fModelPol1 -> GetModelAPosterioriProbability()) 
    fModelPol2 -> DoGoodnessOfFitTest(100, fModelPol2 -> GetBestFitParameters(), grid, limits) -> 
      Print("modelpol2_gof.ps", 2, TMath::Log10(fModelPol2 -> Likelihood(fModelPol2 -> GetBestFitParameters())));

  // ---------------------------------------------------------
  // summarize
  // ---------------------------------------------------------
  
  fModelManager -> PrintSummary(); 

  fModelManager -> PrintSummaryToFile("output.log"); 

  // ---------------------------------------------------------
  // Print data with best fit result 
  // ---------------------------------------------------------

  TCanvas* canvas_bestfit = new TCanvas(); 
  canvas_bestfit -> cd(); 

  TH2D* hist_axes = new TH2D("hist_axes", "", 1, 0.0, 100.0, 1, 0.0, 6.0); 
  hist_axes -> SetXTitle("x"); 
  hist_axes -> SetYTitle("y"); 
  hist_axes -> SetStats(false); 
  hist_axes -> Draw(); 

  TGraphErrors* graph = new TGraphErrors(); 
  graph -> SetMarkerStyle(20); 
  graph -> SetMarkerColor(kBlack); 

  for (int i = 0; i < fModelManager -> GetNDataPoints(); i++)
    {
      graph -> SetPoint(i, 
			fModelManager -> GetDataPoint(i) -> GetValue(0), 
			fModelManager -> GetDataPoint(i) -> GetValue(1)); 
      graph -> SetPointError(i, 
			     0.0, 
			     fModelManager -> GetDataPoint(i) -> GetValue(2)); 
    }
  
  TF1* func_pol0 = new TF1("func_pol0", "[0]", 0, 100); 
  func_pol0 -> SetLineWidth(3); 
  func_pol0 -> SetLineColor(kBlack); 
  
  func_pol0 -> SetParameter(0, fModelPol0 -> GetBestFitParameter(0)); 

  func_pol0 -> Draw("SAME"); 

  TF1* func_pol1 = new TF1("func_pol1", "[0] + [1] * x", 0, 100); 
  func_pol1 -> SetLineWidth(3); 
  func_pol1 -> SetLineColor(kRed); 

  func_pol1 -> SetParameter(0, fModelPol1 -> GetBestFitParameter(0)); 
  func_pol1 -> SetParameter(1, fModelPol1 -> GetBestFitParameter(1)); 

  func_pol1 -> Draw("SAME"); 

  TF1* func_pol2 = new TF1("func_pol2", "[0] + [1] * x + [2] * x * x", 0, 100); 
  func_pol2 -> SetLineWidth(3); 
  func_pol2 -> SetLineColor(kGreen); 

  func_pol2 -> SetParameter(0, fModelPol2 -> GetBestFitParameter(0)); 
  func_pol2 -> SetParameter(1, fModelPol2 -> GetBestFitParameter(1)); 
  func_pol2 -> SetParameter(2, fModelPol2 -> GetBestFitParameter(2)); 

  func_pol2 -> Draw("SAME"); 

  graph -> Draw("SAMEP"); 

  TLegend* legend = new TLegend(0.65, 0.70, 0.95, 0.95); 
  legend -> SetFillColor(kWhite); 
  legend -> SetBorderSize(0); 
  legend -> AddEntry(func_pol0, "ModelPol0", "L"); 
  legend -> AddEntry(func_pol1, "ModelPol1", "L"); 
  legend -> AddEntry(func_pol2, "ModelPol2", "L"); 
  legend -> Draw("SAME"); 

  canvas_bestfit -> Print("data.ps"); 

  // ---------------------------------------------------------
  // close log file 
  // ---------------------------------------------------------

  BCLog::CloseLog(); 

  return 0; 

}

// ---------------------------------------------------------
  
