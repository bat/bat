#include <BCModelPol1.h>
#include <BCLog.h>
#include <BCModelOutput.h> 

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
	// ----------------------------------------------------------

	// calls a function which defines a nicer style than the ROOT
	// default.
	SetStyle(); 

	// ---------------------------------------------------------
	// open log file 
	// ---------------------------------------------------------

	// opens the log file. 
	BCLog::OpenLog("log.txt", BCLog::detail, BCLog::detail); 

	// ---------------------------------------------------------
	// model definition 
	// ---------------------------------------------------------

	// creates a new model of type BCModelPol1 with the name
	// "ModelPol1". the model is defined in the two files
	// src/BCModelPol1.h and src/BCModelPol1.cxx. 
	BCModelPol1* fModelPol1 = new BCModelPol1("ModelPol1"); 

	// add parameters to the model. the first parameter will have the
	// index 0, the seoncd parameter will have the index 1.
	fModelPol1 -> AddParameter("constant",  1.0,  3.0);  // index 0
	fModelPol1 -> AddParameter("slope",    -0.03, 0.03); // index 1

	// ---------------------------------------------------------
	// model output file 
	// ---------------------------------------------------------
	
	// creates a ROOT output file which stores all the necessary
	// information.
	BCModelOutput * fModelOutputPol1 = new BCModelOutput(fModelPol1, "output.root"); 

	// the actual markov chain will be written to the output file. this
	// might result in a large file. 
	fModelOutputPol1 -> WriteMarkovChain(true); 
	
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

	// the allowed range of data values has to be defined for error
	// propagation and fitting. 
	fModelPol1 -> SetDataBoundaries(0, 0.0, 100.0); // possible x-values
	fModelPol1 -> SetDataBoundaries(1, 0.0, 5.0); // possible y-values
	fModelPol1 -> SetDataBoundaries(2, 0.2, 5.2); // possible sigmas

	// in this example the data points are not chosen arbitrarily but
	// they have fixed values on the x-axis, i.e, 5, 15, ... . also, the
	// resolution is fixed. in general this can be solved by fixing the
	// "axes".
	fModelPol1 -> FixDataAxis(0, true); // fixes x-values 
	fModelPol1 -> FixDataAxis(2, true); // fixes resolution 

	// ---------------------------------------------------------
	// optimization 
	// ---------------------------------------------------------

	// set the optimization algorithm to the ROOT-version of Minuit. 
	fModelPol1 -> SetOptimizationMethod(BCIntegrate::kOptMinuit); 

	// perform optimization. the best fit parameter values will be
	// stored and can later be retrieved.
	fModelPol1 -> FindMode(); 

	// ---------------------------------------------------------
	// marginalize 
	// ---------------------------------------------------------

	// during the marginalization, the error propagation is done. thus,
	// the x- and y-indices of the data values have to be set.
	fModelPol1 -> SetFitFunctionIndices(0, 1); 

	// marginalizes the probability density with respect to all
	// parameters, i.e. constant and slope, and with respect to all
	// combinations of two parameters, in this case constant-slope.
	fModelPol1 -> MarginalizeAll();

	// the one- and two-dimensional marginalized probability densities
	// are kept in memory and are returned from the model class. they
	// are printed into a .ps file.
	fModelPol1 -> PrintAllMarginalized("plots.ps"); 

	// ---------------------------------------------------------
	// Do goodness-of-fit test  
	// ---------------------------------------------------------

	// once the most probable model was found the goodness-of-fit can be
	// tested. the question to be answered is: given the best fit
	// parameters what is the probability to obtain a better agreement
	// (conditional probability)? the answer can be given by integrating
	// over data. 

	// calculate the p-value for the set of best-fit parameters.
	//	fModelPol1 -> CalculatePValue(fModelPol1 -> GetBestFitParameters()); 
	fModelPol1 -> CalculatePValue(fModelPol1 -> GetBestFitParameters()); 
	
	// ---------------------------------------------------------
	// write to output file 
	// ---------------------------------------------------------

	// fill the output of the analysis into the ROOT file. 
	fModelOutputPol1 -> FillAnalysisTree(); 

	// fill the Markov Chain into the ROOT file. 
	fModelOutputPol1 -> WriteMarginalizedDistributions(); 

	// write to file and close 
	fModelOutputPol1 -> Close(); 

	// ---------------------------------------------------------
	// summarize
	// ---------------------------------------------------------

	// prints the results of the analysis into a file 
	fModelPol1 -> PrintResults("summary.txt"); 

	// ---------------------------------------------------------
	// Print data with best fit result 
	// ---------------------------------------------------------
	
	// defines a new canvas. 
	TCanvas* canvas_bestfit = new TCanvas(); 
	canvas_bestfit -> cd(); 

	// defines a histogram for the axes and draws it. 
	TH2D* hist_axes = new TH2D("hist_axes", "Data;x;y", 1, 0.0, 100.0, 1, 1.5, 3.5); 
	hist_axes -> SetStats(false); 
	hist_axes -> Draw(); 

	// draw the error band 
	fModelPol1 -> GetErrorBandGraph(0.16, 0.84) -> Draw("SAMEF"); 

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
	fModelPol1 -> GetFitFunctionGraph(fModelPol1 -> GetBestFitParameters()) -> Draw("SAMEC"); 
	
	// draw the graph containing the data. 
	graph -> Draw("SAMEP"); 

	// print the canvas to a .ps file 
	canvas_bestfit -> Print("fit.ps"); 
	
	// ---------------------------------------------------------
	// close log file 
	// ---------------------------------------------------------

	// closes the log file 
	BCLog::CloseLog(); 

	return 1; 

}

// ---------------------------------------------------------
  
