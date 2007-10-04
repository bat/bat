#include <BCModelPol1.h>
#include <BCLog.h>

#include "style.c" 

#include <TCanvas.h> 
#include <TGraphErrors.h>
#include <TH2D.h> 
#include <TF1.h> 

// ---------------------------------------------------------
  
int main()
{

	// ---------------------------------------------------------
	// set style  
	// --------------------- ------------------------------------

	// calls a function which defines a nice style. 

	SetStyle(); 

	// ---------------------------------------------------------
	// open log file 
	// ---------------------------------------------------------

	// opens the log file. 

	BCLog::OpenLog(); 

	// ---------------------------------------------------------
	// model definition 
	// ---------------------------------------------------------

	// creates a new model of type BCModelPol1 with the name
	// "ModelPol1". the model is defined in the two files
	// src/BCModelPol1.h and src/BCModelPol1.cxx. The two parameters of
	// the model and their ranges are defined automatically by
	// construction of the model.

	BCModelPol1* fModelPol1 = new BCModelPol1("ModelPol1"); 

	// ---------------------------------------------------------
	// read data from file 
	// ---------------------------------------------------------

	// creates a new data set. 

	BCDataSet* fDataSet = new BCDataSet(); 

	// data is read in from a text file. three values per data point are
	// read in: x, y, and the uncertainty on y. if the file is not found
	// or corrupt, the program returns -1.

	if (fDataSet -> ReadDataFromFileTxt("./data/data.txt", 3) != 0)
		return -1; 

	// assigns the data set to the model defined previously. 

	fModelPol1 -> SetDataSet(fDataSet); 

	// ---------------------------------------------------------
	// find mode 
	// ---------------------------------------------------------

	// runs the algorithm which searches the whole parameter space for
	// the maximum probability (mode). for details on the algorithm
	// consult the manual. this step might need a while, depending on
	// the function.

	fModelPol1 -> SetModeFindingMethod(BCIntegrate::kMFMinuit); 

	fModelPol1 -> FindMode(); 

	// ---------------------------------------------------------
	// prepare calculation of error band 
	// ---------------------------------------------------------

	// set limits on data values 

	fModelPol1 -> SetDataBoundaries(0, 0.0, 100.0); 
	fModelPol1 -> SetDataBoundaries(1, 1.0,   4.0); 
	fModelPol1 -> SetDataBoundaries(2, 0.2,   0.2); 

	// set x and y value indices 

	fModelPol1 -> SetFitFunctionIndices(0, 1); 

	// ---------------------------------------------------------
	// marginalize 
	// ---------------------------------------------------------

	// marginalizes the probability density with respect to all
	// parameters, i.e. constant and slope and with respect to all
	// combinations of two parameters, in this case constant-slope. the
	// number of bins define the numerical precision. 

	fModelPol1 -> SetNbins(100);
	fModelPol1 -> MarginalizeAll();

	// the one-dimensional marginalized probability densities are kept
	// in memory and are returned from the model class. they are printed
	// into a .ps file. 

	fModelPol1 -> GetMarginalized("constant") -> Print("modelpol1_constant.ps");
	fModelPol1 -> GetMarginalized("slope")    -> Print("modelpol1_slope.ps");

	// the two-dimensional marginalized probability densitiy is kept in
	// memory and is returned from the model class. it is printed into a
	// .ps file.

	fModelPol1 -> GetMarginalized("constant", "slope") -> Print("modelpol1_constant_slope.ps", 2);

	// ---------------------------------------------------------
	// summarize
	// ---------------------------------------------------------

	// prints a summary of the model, the parameter estimate, etc. on
	// the screen and to a file.

	fModelPol1 -> PrintSummary(); 

	// ---------------------------------------------------------
	// Print data with best fit result 
	// ---------------------------------------------------------

	// defines a new canvas. 

	TCanvas* canvas_bestfit = new TCanvas(); 
	canvas_bestfit -> cd(); 

	// defines a histogram for the axes and draws it. 

	TH2D* hist_axes = new TH2D("hist_axes", "Data;x;y", 1, 0.0, 100.0, 1, 0.0, 6.0); 
	hist_axes -> SetStats(false); 
	hist_axes -> Draw(); 

	// draw the error band 

	fModelPol1 -> GetErrorBandXY() -> Draw("COL"); 
	//	fModelPol1 -> GetErrorBandGraph(0.16, 0.84) -> Draw("F"); 

	// defines a graph with errors. 

	TGraphErrors* graph = new TGraphErrors(); 
	graph -> SetMarkerStyle(20); 

	// sets the points of the graph to the data points read in
	// previously. loops over all entries. 

	for (int i = 0; i < fModelPol1 -> GetNDataPoints(); i++)
		{
			graph -> SetPoint(i, 
												fModelPol1 -> GetDataPoint(i) -> GetValue(0), 
												fModelPol1 -> GetDataPoint(i) -> GetValue(1)); 
			graph -> SetPointError(i, 
														 0.0, 
														 fModelPol1 -> GetDataPoint(i) -> GetValue(2)); 
		}

	// draw best fit function 

	fModelPol1 -> GetFitFunctionGraph(fModelPol1 -> GetBestFitParametersMarginalized()) -> Draw("SAMEC"); 

	// draw the graph containing the data. 

	graph -> Draw("SAMEP"); 

	// print the canvas to a .ps file 

	canvas_bestfit -> Print("data.ps"); 

	// ---------------------------------------------------------
	// close log file 
	// ---------------------------------------------------------

	// closes the log file 

	BCLog::CloseLog(); 

	return 0; 

}

// ---------------------------------------------------------
  
