#include <BCModelPol0.h>
#include <BCModelPol1.h>
#include <BCModelPol2.h>
#include <BCModelManager.h>
#include <BCDataPoint.h> 
#include <BCLog.h> 

#include "style.c"

#include <TCanvas.h> 
#include <TGraphErrors.h>
#include <TH2D.h> 
#include <TF1.h> 
#include <TLegend.h> 

#include <iostream>

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

	// creates three new models which describe polynomials of order 0, 1
	// and 2. the are defined in the /src directory. the number of
	// parameters equals the order of the polynomial. they are defined
	// automatically by construction of the particular model.

	BCModelPol0* fModelPol0 = new BCModelPol0("ModelPol0"); 
	BCModelPol1* fModelPol1 = new BCModelPol1("ModelPol1"); 
	BCModelPol2* fModelPol2 = new BCModelPol2("ModelPol2"); 

	// ---------------------------------------------------------
	// define model manager 
	// ---------------------------------------------------------

	// defines a model manager. this manager will be used to make sure
	// all models use the same data set. 

	BCModelManager* fModelManager = new BCModelManager(); 

	// adds all three models to the manager and passes the a priori
	// probability to the manager. 

	fModelManager -> AddModel(fModelPol0, 1.0/3.0); 
	fModelManager -> AddModel(fModelPol1, 1.0/3.0); 
	fModelManager -> AddModel(fModelPol2, 1.0/3.0); 

	// ---------------------------------------------------------
	// read data from file 
	// ---------------------------------------------------------

	// creates a new data set 

	BCDataSet* fDataSet = new BCDataSet(); 

	// data is read in from a text file. three values per data point are
	// read in: x, y, and the uncertainty on y. if the file is not found
	// or corrupt, the program returns -1.

	if (fDataSet -> ReadDataFromFileTxt("./data/data_ModelPol2.txt", 3) != 0)
		return -1;

	// assigns the data set to the model manager. the manager
	// automatically assigns the same data set to all three models. 

	fModelManager -> SetDataSet(fDataSet); 

	// ---------------------------------------------------------
	// prepare calculation of error band and goodness-of-fit test 
	// ---------------------------------------------------------

	// sets boundaries on possible data values. this is needed since the
	// integration is performed over the data values. in this example
	// the uncertainty on the measured y value is fixed.

	// set limits on data values 

	fModelPol0 -> SetDataBoundaries(0, 0.0, 100.0); 
	fModelPol0 -> SetDataBoundaries(1, 0.0,   4.5); 
	fModelPol0 -> SetDataBoundaries(2, 0.2,   0.2); 

	fModelPol1 -> SetDataBoundaries(0, 0.0, 100.0); 
	fModelPol1 -> SetDataBoundaries(1, 0.0,   4.5); 
	fModelPol1 -> SetDataBoundaries(2, 0.2,   0.2); 

	fModelPol2 -> SetDataBoundaries(0, 0.0, 100.0); 
	fModelPol2 -> SetDataBoundaries(1, 0.0,   4.5); 
	fModelPol2 -> SetDataBoundaries(2, 0.2,   0.2); 

	// set x and y value indices 

	fModelPol0 -> SetFitFunctionIndices(0, 1); 
	fModelPol1 -> SetFitFunctionIndices(0, 1); 
	fModelPol2 -> SetFitFunctionIndices(0, 1); 

	// ---------------------------------------------------------
	// initialize 
	// ---------------------------------------------------------

	// initializes the manager. it automatically calculates the
	// normalizations for each model and calculates the model a
	// posteriori probabilities. this step might take a while, depending
	// on the number of parameters.


	fModelManager -> Initialize(); 

	fModelPol0 -> SetModeFindingMethod(BCIntegrate::kMFMinuit); 
	fModelPol1 -> SetModeFindingMethod(BCIntegrate::kMFMinuit); 
	fModelPol2 -> SetModeFindingMethod(BCIntegrate::kMFMinuit); 

	fModelPol0 -> FindMode(); 
	fModelPol1 -> FindMode(); 
	fModelPol2 -> FindMode(); 

/* 	// --------------------------------------------------------- */
/* 	// marginalize  */
/* 	// --------------------------------------------------------- */

	// selects the most probable model and marginalizes the probability
	// density with respect to all parameters

	// first, remember the model a posterior probabilities. 

	double post_modelpol0 = fModelPol0 -> GetModelAPosterioriProbability(); 
	double post_modelpol1 = fModelPol1 -> GetModelAPosterioriProbability(); 
	double post_modelpol2 = fModelPol2 -> GetModelAPosterioriProbability(); 

	// select 0th order polynomial model. 

	// marginalizes the probability density with respect to all
	// parameters, i.e. constant and slope and with respect to all
	// combinations of two parameters, in this case constant-slope. the
	// number of bins define the numerical precision. 
	
	fModelPol0 -> SetNbins(100); 
	fModelPol0 -> MarginalizeAll();
	
	// the one-dimensional marginalized probability densities are kept
	// in memory and are returned from the model class. they are printed
	// into a .ps file. 
	
	fModelPol0 -> GetMarginalized("constant") -> Print("modelpol0_constant.ps");

	// select 1st order polynomial model. 
	
	// marginalizes the probability density with respect to all
	// parameters, i.e. constant and slope and with respect to all
	// combinations of two parameters, in this case constant-slope. the
	// number of bins define the numerical precision. 
	
	fModelPol1 -> SetNbins(100); 
	fModelPol1 -> MarginalizeAll();
	
	// the one-dimensional marginalized probability densities are kept
	// in memory and are returned from the model class. they are printed
	// into a .ps file. 
	
	fModelPol1 -> GetMarginalized("constant")          -> Print("modelpol1_constant.ps");
	fModelPol1 -> GetMarginalized("slope")             -> Print("modelpol1_slope.ps");
	
	// the two-dimensional marginalized probability densitiy is kept in
	// memory and is returned from the model class. it is printed into a
	// .ps file.
	
	fModelPol1 -> GetMarginalized("constant", "slope") -> Print("modelpol1_constant_slope.ps", 2);

	// select 2nd order polynomial model. 

	// marginalizes the probability density with respect to all
	// parameters, i.e. constant and slope and with respect to all
	// combinations of two parameters, in this case constant-slope. the
	// number of bins define the numerical precision. 
	
	fModelPol2 -> SetNbins(100); 
	fModelPol2 -> MarginalizeAll();
	
	// the one-dimensional marginalized probability densities are kept
	// in memory and are returned from the model class. they are printed
	// into a .ps file. 
	
	fModelPol2 -> GetMarginalized("constant")          -> Print("modelpol2_constant.ps");
	fModelPol2 -> GetMarginalized("slope")             -> Print("modelpol2_slope.ps");
	fModelPol2 -> GetMarginalized("quad")              -> Print("modelpol2_quad.ps");
	
	// the one-dimensional marginalized probability densities are kept
	// in memory and are returned from the model class. they are printed
	// into a .ps file. 
	
	fModelPol2 -> GetMarginalized("constant", "slope") -> Print("modelpol2_constant_slope.ps", 2);
	fModelPol2 -> GetMarginalized("constant", "quad")  -> Print("modelpol2_constant_quad.ps", 2);
	fModelPol2 -> GetMarginalized("slope",    "quad")  -> Print("modelpol2_slope_quad2.ps", 2);

	// ---------------------------------------------------------
	// Do goodness-of-fit test  
	// ---------------------------------------------------------

	// once the most probable model was found the goodness-of-fit can be
	// tested. the question to be answered is: given the best fit
	// parameters what is the probability to obtain a better agreement
	// (conditional probability)? the answer can be given by integrating
	// over data. technically this is done by creating ensembles given
	// the best fit parameters. 

	// in this example the data points are not chosen arbitrarily but
	// they have defined values on a grid, e.g. 5, 15, ... in general
	// this can be solved by defining a grid for the values which are
	// not free.

	// defines a vector which contains the information if a value is
	// defined on a grid (true) or is free (false). in this example only
	// the x component is defined on a grid, y and the uncertainty on y
	// are free (keep in mind that the uncertainty on y is a fixed
	// number, see previously defined boundaries). 

	std::vector <bool> grid;
	grid.push_back(true);
	grid.push_back(false);
	grid.push_back(false);

	// defines a vector which describes the grid. the first entry
	// defines the first value, e.g. 5, the second entry defines the
	// step size, e.g. 10, the third value defines the number of grid
	// points, e.g. 10. if more than one variable is defined on a grid
	// the same sequence is push back into the vector for the other
	// variables. 

	std::vector <double> limits;
	limits.push_back( 5.0);
	limits.push_back(10.0);
	limits.push_back(10.0);

	// selects the most probable model and performs a goodness-of-fit
	// test. ensembles of data sets are created, given the best fit
	// parameters, and their conditional probability is evaluated. the
	// values are histogrammed and the frequency distribution is
	// interpreted as a probability density for the agreement. the
	// p-value is estimated. for details on the goodness-of-fit test,
	// see the manual. 

	//	if (post_modelpol0 > post_modelpol1 && 
	//			post_modelpol0 > post_modelpol2) 
	//		{
			// performs the goodness-of-fit test. creates 100 ensembles
			// given the best fit parameters. the frequency distribution is
			// printed into a .ps file and the conditional probability for
			// the original data is indicated by a line. 

			fModelPol0 -> DoGoodnessOfFitTest(200, fModelPol0 -> GetBestFitParameters(), grid, limits) -> 
				Print("modelpol0_gof.ps", 1, TMath::Log10(fModelPol0 -> Likelihood(fModelPol0 -> GetBestFitParameters())));
			//		}

			//	if (post_modelpol1 > post_modelpol0 && 
			//			post_modelpol1 > post_modelpol2) 
			//		{
			// performs the goodness-of-fit test. creates 100 ensembles
			// given the best fit parameters. the frequency distribution is
			// printed into a .ps file and the conditional probability for
			// the original data is indicated by a line. 

			fModelPol1 -> DoGoodnessOfFitTest(200, fModelPol1 -> GetBestFitParameters(), grid, limits) -> 
				Print("modelpol1_gof.ps", 1, TMath::Log10(fModelPol1 -> Likelihood(fModelPol1 -> GetBestFitParameters())));
			//		}

			//	if (post_modelpol2 > post_modelpol0 && 
			//			post_modelpol2 > post_modelpol1) 
			//		{
			// performs the goodness-of-fit test. creates 100 ensembles
			// given the best fit parameters. the frequency distribution is
			// printed into a .ps file and the conditional probability for
			// the original data is indicated by a line. 

			fModelPol2 -> DoGoodnessOfFitTest(200, fModelPol2 -> GetBestFitParameters(), grid, limits) -> 
				Print("modelpol2_gof.ps", 1, TMath::Log10(fModelPol2 -> Likelihood(fModelPol2 -> GetBestFitParameters())));
			//		}

	// ---------------------------------------------------------
	// summarize
	// ---------------------------------------------------------

	// prints a summary of the model, the parameter estimate, etc. on
	// the screen and to a file.

	fModelManager -> PrintSummary(); 

	fModelManager -> PrintSummary("output.log");

	// ---------------------------------------------------------
	// Print data with best fit result 
	// ---------------------------------------------------------

	// defines a new canvas 

	TCanvas* canvas_summary = new TCanvas(); 
	canvas_summary -> cd(); 

	// defines a graph with errors 

	TGraphErrors* graph_data = new TGraphErrors(); 
	graph_data -> SetMarkerStyle(20); 
	graph_data -> SetMarkerColor(kBlack); 

	// sets the points of the graph to the data points read in
	// previously. loops over all entries. 

	for (int i = 0; i < fModelManager -> GetNDataPoints(); i++)
		{
			graph_data -> SetPoint(i, 
														 fModelManager -> GetDataPoint(i) -> GetValue(0), 
														 fModelManager -> GetDataPoint(i) -> GetValue(1)); 
			graph_data -> SetPointError(i, 
																	0.0, 
																	fModelManager -> GetDataPoint(i) -> GetValue(2)); 
		}

	// defines a histogram for the axes and draws it.

	TH2D* hist_axes = new TH2D("hist_axes", "", 1, 0.0, 100.0, 1, 0.0, 6.0);
	hist_axes -> SetXTitle("x");
	hist_axes -> SetYTitle("y");
	hist_axes -> SetStats(false);
	hist_axes -> Draw();

	graph_data -> Draw("SAMEP"); 

	// define best fit function graphs 

	TGraph * graph_bestfit_pol0 = fModelPol0 -> GetFitFunctionGraph(); 
	graph_bestfit_pol0 -> SetLineColor(kRed); 

	TGraph * graph_bestfit_pol1 = fModelPol1 -> GetFitFunctionGraph(); 
	graph_bestfit_pol1 -> SetLineColor(kBlack); 

	TGraph * graph_bestfit_pol2 = fModelPol2 -> GetFitFunctionGraph(); 
	graph_bestfit_pol2 -> SetLineColor(kGreen); 

	graph_bestfit_pol0 -> Draw("SAMEC"); 
	graph_bestfit_pol1 -> Draw("SAMEC"); 
	graph_bestfit_pol2 -> Draw("SAMEC"); 
			
	// define and draw a legend. 

	TLegend* legend = new TLegend(0.65, 0.70, 0.95, 0.95); 
	legend -> SetFillColor(kWhite); 
	legend -> SetBorderSize(0); 
	legend -> AddEntry(graph_data, "Data", "P"); 
	legend -> AddEntry(graph_bestfit_pol0, "ModelPol0", "L"); 
	legend -> AddEntry(graph_bestfit_pol1, "ModelPol1", "L"); 
	legend -> AddEntry(graph_bestfit_pol2, "ModelPol2", "L"); 
	legend -> Draw("SAME"); 

	// print the canvas to a .ps file 

	canvas_summary -> Print("data_allmodels.ps"); 

	// defines a new canvas 

	TCanvas* canvas_bestfit = new TCanvas(); 
	canvas_bestfit -> cd(); 

	// draw error band and best fit function of best model 

	if (post_modelpol0 > post_modelpol1 && 
			post_modelpol0 > post_modelpol2) 
		{
			hist_axes -> Draw();
			fModelPol0 -> GetErrorBandGraph(0.16, 0.84) -> Draw("F"); 
			//			fModelPol0 -> GetErrorBandXY() -> Draw("COL");	
			graph_bestfit_pol0 -> SetLineColor(kBlack); 
			graph_bestfit_pol0 -> Draw("SAMEC"); 
		}

	if (post_modelpol1 > post_modelpol0 && 
			post_modelpol1 > post_modelpol2) 
		{
			hist_axes -> Draw();
			fModelPol1 -> GetErrorBandGraph(0.16, 0.84) -> Draw("F"); 
			//			fModelPol1 -> GetErrorBandXY() -> Draw("COL");	
			graph_bestfit_pol1 -> SetLineColor(kBlack); 
			graph_bestfit_pol1 -> Draw("SAMEC"); 
		}

	if (post_modelpol2 > post_modelpol0 && 
			post_modelpol2 > post_modelpol1) 
		{
			hist_axes -> Draw();
			fModelPol2 -> GetErrorBandGraph(0.16, 0.84) -> Draw("F"); 
			//			fModelPol2 -> GetErrorBandXY() -> Draw("COL");	
			graph_bestfit_pol2 -> SetLineColor(kBlack); 
			graph_bestfit_pol2 -> Draw("SAMEC"); 
		}

	// draw data points 

	graph_data -> Draw("SAMEP"); 

	// print the canvas to a .ps file 

	canvas_bestfit -> Print("data_bestfit.ps"); 

	// ---------------------------------------------------------
	// close log file 
	// ---------------------------------------------------------

	// closes the log file. 

	BCLog::CloseLog(); 

	return 0; 

}

// ---------------------------------------------------------
  
