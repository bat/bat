#include <BCModelEfficiency.h>
#include <BCLog.h>
#include <BCModelOutput.h> 
#include <BCDataSet.h>
#include <BCAux.h>

#include <TCanvas.h> 
#include <TGraphErrors.h>
#include <TH2D.h> 
#include <TF1.h> 

// ---------------------------------------------------------
// This examples calculates an efficiency based on two numbers and the
// binomial distribution. 
// ---------------------------------------------------------
  
int main()
{
  
	// ---------------------------------------------------------
	// set style  
	// ----------------------------------------------------------

	// calls a function which defines a nicer style than the ROOT default
	BCAux::SetStyle();

	// ---------------------------------------------------------
	// open log file 
	// ---------------------------------------------------------

	// opens the log file. All details of the calculation are printed to
	// the screen and written to the log file.
	BCLog::OpenLog("log.txt", BCLog::detail, BCLog::detail); 

	// ---------------------------------------------------------
	// model definition 
	// ---------------------------------------------------------

	// creates a new model of type BCModelEfficiency with the name
	// "ModelEfficiency". the model is defined in the two files
	// src/BCModelEfficiency.h and src/BCModelEfficiency.cxx.
	BCModelEfficiency* fModelEfficiency = new BCModelEfficiency("ModelEfficiency"); 

	// adds the parameter "epsilon" to the model. 
	fModelEfficiency -> AddParameter("epsilon",  0.0,  1.0);  // index 0

	// ---------------------------------------------------------
	// model output file 
	// ---------------------------------------------------------
	
	// creates a ROOT output file which stores all the necessary
	// information.
	BCModelOutput * fModelOutputEfficiency = new BCModelOutput(fModelEfficiency, "output.root"); 

	// the actual markov chain will be written to the output file. this
	// might result in a large file. 
	//	fModelOutputEfficiency -> WriteMarkovChain(true); 
	
	// ---------------------------------------------------------
	// read data from file 
	// ---------------------------------------------------------

	// creates a new data set. 
	BCDataSet* fDataSet = new BCDataSet(); 

	// data is read in from a text file. 
	if (fDataSet -> ReadDataFromFileTxt("./data/data.txt", 2) != 0)
		return -1; 

	// assigns the data set to the model defined previously.
	fModelEfficiency -> SetDataSet(fDataSet); 

	// the allowed range of data values has to be defined for error
	// propagation and fitting. 
	double Nmax = fDataSet -> GetDataPoint(0) -> GetValue(1); 
	fModelEfficiency -> SetDataBoundaries(0, 0.0, Nmax, false);
	fModelEfficiency -> SetDataBoundaries(1, Nmax, Nmax); 

	// Fix the maximum number. 
	fModelEfficiency -> FixDataAxis(1, true); 

	// ---------------------------------------------------------
	// optimization 
	// ---------------------------------------------------------

	// set the optimization algorithm to the ROOT-version of Minuit. 
	fModelEfficiency -> SetOptimizationMethod(BCIntegrate::kOptMinuit); 

	// perform optimization. the best fit parameter values will be
	// stored and can later be retrieved.
	fModelEfficiency -> FindMode(); 

	// ---------------------------------------------------------
	// marginalize 
	// ---------------------------------------------------------

	// marginalizes the probability density. in this case, the
	// distribution for the parameter "epsilon" will be created. 
	fModelEfficiency -> MarginalizeAll();

	// the (marginalized) probability density are kept in memory and
	// are returned from the model class. they are printed into a .ps
	// file.
	fModelEfficiency -> PrintAllMarginalized("plots.ps"); 

	// ---------------------------------------------------------
	// Do goodness-of-fit test  
	// ---------------------------------------------------------

	// once the most probable model was found the goodness-of-fit can be
	// tested. the question to be answered is: given the best fit
	// parameters what is the probability to obtain a better agreement
	// (conditional probability)? the answer can be given by integrating
	// over data. 

	// calculate the p-value for the set of best-fit parameters.
	fModelEfficiency -> CalculatePValue(fModelEfficiency -> GetBestFitParameters()); 
	
	// ---------------------------------------------------------
	// write to output file 
	// ---------------------------------------------------------

	// fill the output of the analysis into the ROOT file. 
	fModelOutputEfficiency -> FillAnalysisTree(); 

	// fill the Markov Chain into the ROOT file. 
	fModelOutputEfficiency -> WriteMarginalizedDistributions(); 

	// write to file and close 
	fModelOutputEfficiency -> Close(); 

	// ---------------------------------------------------------
	// summarize
	// ---------------------------------------------------------

	// prints the results of the analysis into a file 
	fModelEfficiency -> PrintResults("summary.txt"); 
	
	// ---------------------------------------------------------
	// close log file 
	// ---------------------------------------------------------

	// closes the log file 
	BCLog::CloseLog(); 

	return 1; 

}

// ---------------------------------------------------------
  
