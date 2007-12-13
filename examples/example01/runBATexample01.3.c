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
	// model output file 
	// ---------------------------------------------------------
	
	// creates a ROOT output file which stores all the necessary
	// information. 

	BCModelOutput * fModelOutputPol1 = new BCModelOutput(fModelPol1, "output.root"); 
	
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
	fModelPol1 -> GetMarginalized("constant", "slope") -> Print("modelpol1_constant_slope_color.ps", 1);
	
	// ---------------------------------------------------------
	// write to output file 
	// ---------------------------------------------------------

	// fill the ROOT file with the actual output of the model. 

	fModelOutputPol1 -> FillAnalysisTree(); 

	fModelOutputPol1 -> WriteMarginalizedDistributions(); 

	// write to file and close 

	fModelOutputPol1 -> Close(); 

	// ---------------------------------------------------------
	// summarize
	// ---------------------------------------------------------

	// prints a summary of the model, the parameter estimate, etc. on
	// the screen and to a file.

	fModelPol1 -> PrintSummary(); 

	// ---------------------------------------------------------
	// close log file 
	// ---------------------------------------------------------

	// closes the log file 

	BCLog::CloseLog(); 

	return 0; 

}

// ---------------------------------------------------------
  
