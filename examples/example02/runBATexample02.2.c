#include <BCModelPol0.h>
#include <BCModelPol1.h>
#include <BCModelPol2.h>
#include <BCModelManager.h>
#include <BCModelOutput.h>
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
	// normalize 
	// ---------------------------------------------------------

	// automatically calculates the normalizations for each model and
	// calculates the model a posteriori probabilities. this step might
	// take a while, depending on the number of parameters.

	fModelManager -> SetIntegrationMethod(BCIntegrate::kICuba); 

	fModelManager -> Normalize(); 

	// ---------------------------------------------------------
	// find mode 
	// ---------------------------------------------------------

	// finds mode for all models 

	fModelManager -> SetModeFindingMethod(BCIntegrate::kMFMinuit); 

	fModelManager -> FindMode(); 

	// ---------------------------------------------------------
	// summarize
	// ---------------------------------------------------------

	// prints a summary of the model, the parameter estimate, etc. on
	// the screen and to a file.

	fModelManager -> PrintSummary(); 

	fModelManager -> PrintSummary("output.log");

	// ---------------------------------------------------------
	// close log file 
	// ---------------------------------------------------------

	// closes the log file. 

	BCLog::CloseLog(); 

	return 0; 

}

// ---------------------------------------------------------
  
