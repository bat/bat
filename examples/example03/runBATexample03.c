#include <BCModelBackground.h>
#include <BCModelSignal.h>
#include <BCParameter.h>
#include <BCModelManager.h>
#include <BCModelOutput.h>
#include <BCDataPoint.h> 
#include <BCH1D.h> 

#include <TCanvas.h> 
#include <TH1D.h> 
#include <TF1.h> 
#include <BCLog.h> 
#include <TLegend.h> 

#include  <iostream>

#include "style.c"

// ---------------------------------------------------------
  
int main()
{

	// ---------------------------------------------------------
	// set style  
	// ---------------------------------------------------------

	// calls a function which defines a nice style. 

	SetStyle(); 

	// ---------------------------------------------------------
	// open log file 
	// ---------------------------------------------------------

	// opens the log file. 

	BCLog::OpenLog(); 

	// ---------------------------------------------------------
	// define models 
	// ---------------------------------------------------------

	// creates the two models to be tested. The first model describes a
	// flat (background) spectrum while the second one describes a
	// spectrum with a flat (background) and Gaussian (signal)
	// contribution. 

	BCModelBackground* fModelBackground = new BCModelBackground("background model"); 
	BCModelSignal* fModelSignal = new BCModelSignal("signal model"); 

	fModelSignal -> MCMCSetFlagPCA(true); 
	fModelBackground -> MCMCSetTrialFunctionScale(0.01); 
	//	fModelSignal -> MCMCSetTrialFunctionScale(0.01); 

	// ---------------------------------------------------------
	// define manager 
	// ---------------------------------------------------------

	// defines a model manager. this manager will be used to make sure
	// all models use the same data set. 

	BCModelManager* fModelManager = new BCModelManager(); 

	// adds all three models to the manager and passes the a priori
	// probability to the manager. 

	fModelManager -> AddModel(fModelBackground, 0.5); 
	fModelManager -> AddModel(fModelSignal, 0.5); 

	fModelManager -> SetNChains(10); 

	// ---------------------------------------------------------
	// model output file 
	// ---------------------------------------------------------
	
	// creates a ROOT output file which stores all the necessary
	// information. 

	BCModelOutput * fModelOutputBackground = new BCModelOutput(fModelBackground, "output_Background.root"); 
	fModelOutputBackground -> WriteMarkovChain(true); 

	BCModelOutput * fModelOutputSignal = new BCModelOutput(fModelSignal, "output_Signal.root"); 
	fModelOutputSignal -> WriteMarkovChain(true); 

	// ---------------------------------------------------------
	// read data from file 
	// ---------------------------------------------------------

	// creates a new data set 

	BCDataSet* fDataSet = new BCDataSet(); 

	// read data from a text file. two values per data point are
	// defined, i.e. the energy of the bin and the number of events
	// observed in that bin. the first two entries define the region of
	// interest in the spectrum. 

	if (fDataSet -> ReadDataFromFileTxt("./data/data.txt", 2) != 0) 
		return -1; 

	// assigns the data set to the model manager. the manager
	// automatically assigns the same data set to all three models.

	fModelManager -> SetDataSet(fDataSet);

	// ---------------------------------------------------------
	// prepare calculation of error band 
	// ---------------------------------------------------------

	// set limits on data values 

	fModelManager -> SetDataBoundaries(0, fDataSet -> GetDataPoint(0) -> GetValue(0), fDataSet -> GetDataPoint(0) -> GetValue(1)); 
	fModelManager -> SetDataBoundaries(1, 0.0, 20); 

	// set x and y value indices 

	fModelManager -> SetFitFunctionIndices(0, 1); 

	// ---------------------------------------------------------
	// initialize
	// ---------------------------------------------------------

	// initializes the manager. it automatically calculates the
	// normalizations for each model and calculates the model a
	// posteriori probabilities. this step might take a while, depending
	// on the number of parameters.

	//	fModelManager -> SetIntegrationMethod(BCIntegrate::kIntCuba); 

	//	fModelManager -> Normalize();

	// ---------------------------------------------------------
	// find mode 
	// ---------------------------------------------------------

	// finds mode for all models 

	fModelManager -> SetOptimizationMethod(BCIntegrate::kOptMinuit); 

	fModelManager -> FindMode(); 

	// ---------------------------------------------------------
	// marginalize 
	// ---------------------------------------------------------

	// selects the most probable model and marginalizes the probability
	// density with respect to all parameters

	// first, remember the model a posterior probabilities. 

	//	double post_modelbkg = fModelBackground -> GetModelAPosterioriProbability(); 
	//	double post_modelsgn = fModelSignal -> GetModelAPosterioriProbability(); 

	// marginalize background model 

	// marginalize with respect to all (one) parameter.  the
	// number of bins define the numerical precision. 

	fModelBackground -> SetNbins(100); 
	fModelBackground -> MarginalizeAll();

	// the one-dimensional marginalized probability densities are kept
	// in memory and are returned from the model class. they are printed
	// into a .ps file. 

	fModelBackground -> GetMarginalized("background") -> Print("modelbackground_background.ps");

	// marginalize signal model 

	// marginalize with respect to all parameters.  the
	// number of bins define the numerical precision. 

	fModelSignal -> SetNbins(100); 
	fModelSignal -> MarginalizeAll();

	// the one-dimensional marginalized probability densities are kept
	// in memory and are returned from the model class. they are printed
	// into a .ps file. 

	//	fModelSignal -> GetMarginalized("background") -> Print("modelsignal_background.ps");
	//	fModelSignal -> GetMarginalized("signal") -> Print("modelsignal_signal.ps");
			
	// the two-dimensional marginalized probability densitiy is kept in
	// memory and is returned from the model class. it is printed into a
	// .ps file.

	//	fModelSignal -> GetMarginalized("background", "signal") -> Print("modelsignal_background_signal_contour.ps", 2);

	//	fModelSignal -> GetMarginalized("background", "signal") -> Print("modelsignal_background_signal_color.ps", 5);

	double pvaluesignal = 	fModelSignal -> CalculatePValue(fModelSignal -> GetBestFitParameters()); 

	std::cout << pvaluesignal << std::endl; 

	// ---------------------------------------------------------
	// write to output file 
	// ---------------------------------------------------------

	// fill the ROOT file with the actual output of the model. 

	//	fModelOutputBackground -> FillAnalysisTree(); 
	//	fModelOutputBackground -> WriteMarginalizedDistributions(); 

	//	fModelOutputSignal -> FillAnalysisTree(); 
	//	fModelOutputSignal -> WriteMarginalizedDistributions(); 

	// write to file and close 

	//	fModelOutputBackground -> Close(); 
	//	fModelOutputSignal -> Close(); 

 // ---------------------------------------------------------
 // summarize
 // ---------------------------------------------------------

	// fModelManager -> PrintSummary(); 

	// fModelManager -> PrintSummary("output.log"); 

 // ---------------------------------------------------------
 // Print data with best fit result 
 // ---------------------------------------------------------

	// defines a new canvas 
	/*
	TCanvas* canvas_summary = new TCanvas(); 
	canvas_summary -> cd(); 

	// defines the spectrum histogram 

	int Nbins = int(fModelManager -> GetNDataPoints() - 1); 
	double Emin  = fModelManager -> GetDataPoint(0) -> GetValue(0); 
	double Emax  = fModelManager -> GetDataPoint(0) -> GetValue(1); 
	double dE    = (Emax - Emin) / double(Nbins); 

	TH1D* hist_data = new TH1D("hist_data", "", Nbins, Emin, Emax); 
	hist_data -> SetXTitle("E"); 
	hist_data -> SetYTitle("N"); 
	hist_data -> SetStats(false); 

	// sets the points of the graph to the data points read in
	// previously. loops over all entries. 

	for (int i = 1; i < fModelManager -> GetNDataPoints(); i++)
		hist_data -> Fill(fModelManager -> GetDataPoint(i) -> GetValue(0), 
											fModelManager -> GetDataPoint(i) -> GetValue(1)); 

	// draw the spectrum 

	hist_data -> Draw(); 

	// define a constant function with the best fit parameters of the
	// background model.

	TF1* func_bkg = new TF1("func_bkg", "[0]", Emin, Emax); 
	func_bkg -> SetLineWidth(3); 
	func_bkg -> SetLineColor(kRed); 

	func_bkg -> SetParameter(0, fModelBackground -> GetBestFitParameter(0) / (double(Nbins))); 

	// define a constant function with the best fit parameters of the
	// background model.

	TF1* func_sgn = new TF1("func_sgn", 
													"[0] + [1]/(sqrt(2*3.1416) * [2])*exp(-0.5*((x-[3])/[2])**2)", Emin, Emax); 
	func_sgn -> SetLineWidth(3); 
	func_sgn -> SetLineColor(kGreen); 

	func_sgn -> SetParameter(0, fModelSignal -> GetBestFitParameter(0) / double(Nbins)); 
	func_sgn -> SetParameter(1, fModelSignal -> GetBestFitParameter(1) * dE); 
	func_sgn -> SetParameter(2, 2.0); 
	func_sgn -> SetParameter(3, 2039.0); 

	// draws the functions 

	func_bkg -> Draw("SAME"); 
	func_sgn -> Draw("SAME"); 

	// define and draw a legend. 

	TLegend* legend = new TLegend(0.65, 0.70, 0.90, 0.90); 
	legend -> SetFillColor(kWhite); 
	legend -> SetBorderSize(0); 
	legend -> AddEntry(hist_data, "Data", "H"); 
	legend -> AddEntry(func_bkg, "Background", "L"); 
	legend -> AddEntry(func_sgn, "Signal", "L"); 
	legend -> Draw("SAME"); 

	// print the canvas to a .ps file 

	canvas_summary -> Print("data.ps"); 

	// defines a new canvas 

	TCanvas* canvas_bestfit = new TCanvas(); 
	canvas_bestfit -> cd(); 

	// draw data 

	hist_data -> Draw(); 

	// draw error band and best fit function of best model 

	if (post_modelbkg > post_modelsgn)
		{
			hist_data -> Draw();
			fModelBackground -> GetErrorBandGraph(0.16, 0.84) -> Draw("F"); 
			func_bkg -> SetLineColor(kBlack); 
			func_bkg -> Draw("SAME"); 
			//			fModelPol0 -> GetErrorBandXY() -> Draw("COL");	
			//			graph_bestfit_pol0 -> SetLineColor(kBlack); 
			//			graph_bestfit_pol0 -> Draw("SAMEC"); 
		}

	if (post_modelsgn > post_modelbkg)
		{
			hist_data -> Draw();
			fModelSignal -> GetErrorBandGraph(0.16, 0.84) -> Draw("F"); 
			func_sgn -> SetLineColor(kBlack); 
			func_sgn -> Draw("SAME"); 
			//			fModelPol0 -> GetErrorBandXY() -> Draw("COL");	
			//			graph_bestfit_pol0 -> SetLineColor(kBlack); 
			//			graph_bestfit_pol0 -> Draw("SAMEC"); 
		}

	// draw data 

	hist_data -> Draw("SAME"); 

	// print the canvas to a .ps file 

	canvas_bestfit -> Print("data_bestfit.ps"); 
	*/
 // close log file 

 BCLog::CloseLog(); 

 return 0; 

}

// ---------------------------------------------------------
  
