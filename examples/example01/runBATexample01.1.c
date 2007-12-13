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
  
