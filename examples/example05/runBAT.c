#include <BCModelParticleDecay.h>
#include <BCLog.h>
#include <BCModelOutput.h> 
#include <BCModelManager.h>

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

	BCModelParticleDecay* fModelParticleDecay = new BCModelParticleDecay("ModelParticleDecay"); 

	// ---------------------------------------------------------
	// model output file 
	// ---------------------------------------------------------
	
	// creates a ROOT output file which stores all the necessary
	// information. 

	BCModelOutput * fModelOutputParticleDecay = new BCModelOutput(fModelParticleDecay, "output.root"); 

	// ---------------------------------------------------------
	// read data from file 
	// ---------------------------------------------------------

	// creates a new data set. 

	BCDataSet* fDataSet = new BCDataSet(); 

	// data is read in from a text file. three values per data point are
	// read in: x, y, and the uncertainty on y. if the file is not found
	// or corrupt, the program returns -1.

	if (fDataSet -> ReadDataFromFileTxt("./data/data.txt", 2) != 0)
		return -1; 


	// loop over events 

	int nevents = fDataSet -> GetNDataPoints(); 

	// debug 
	//	nevents = 2; 

	for (int ievent = 0; ievent < nevents; ++ievent)
		{

			// assigns the current event as data set 

			fModelParticleDecay -> SetSingleDataPoint(fDataSet, ievent); 

			// ---------------------------------------------------------
			// find mode 
			// ---------------------------------------------------------

			// runs the algorithm which searches the whole parameter space for
			// the maximum probability (mode). for details on the algorithm
			// consult the manual. this step might need a while, depending on
			// the function.
			
			fModelParticleDecay -> SetModeFindingMethod(BCIntegrate::kMFMinuit); 
			
			fModelParticleDecay -> FindMode(); 

			// ---------------------------------------------------------
			// marginalize 
			// ---------------------------------------------------------

			// marginalizes the probability density with respect to all
			// parameters, i.e. constant and slope and with respect to all
			// combinations of two parameters, in this case constant-slope. the
			// number of bins define the numerical precision. 
			
			fModelParticleDecay -> SetNbins(100);
			fModelParticleDecay -> MarginalizeAll();

			// ---------------------------------------------------------
			// write to output file 
			// ---------------------------------------------------------
			
			// fill the ROOT file with the actual output of the model. 
			
			fModelOutputParticleDecay -> Fill(); 

		}			

	// ---------------------------------------------------------
	// close output file 
	// ---------------------------------------------------------

	// write to file and close 

	fModelOutputParticleDecay -> Close(); 

	// ---------------------------------------------------------
	// close log file 
	// ---------------------------------------------------------

	// closes the log file 

	BCLog::CloseLog(); 

	return 0; 

}

// ---------------------------------------------------------
  
