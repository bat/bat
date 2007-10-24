#include "BCModelOutput.h" 

#include "TDirectory.h"

// --------------------------------------------------------- 

BCModelOutput::BCModelOutput() 
{

	fIndex = 0; 

	fOutputFile = 0; 
	fOutputTree = 0; 
	fModel = 0; 
	fOutputFile = 0; 

	this -> InitializeTree(); 

}

// --------------------------------------------------------- 

BCModelOutput::BCModelOutput(BCModel * model, const char * filename) 
{

	fIndex = 0; 

	fModel = model; 
	
	fFilename = filename; 

	// remeber current directory 

	TDirectory * dir = gDirectory; 

	// create a new file 

	fOutputFile = new TFile(fFilename, "RECREATE"); 

	// change again to the old directory 

	gDirectory = dir; 

	this -> InitializeTree(); 

}

// --------------------------------------------------------- 

BCModelOutput::~BCModelOutput() 
{

	if (fOutputFile)
		{
			fOutputFile -> Close(); 

			delete fOutputFile; 
		} 

	if (fOutputTree)
		delete fOutputTree; 

	if (fOutputFile)
		delete fOutputFile; 

	fModel = 0; 

}

// --------------------------------------------------------- 

void BCModelOutput::SetFile(const char * filename) 
{

	// delete the old file 

	if (fOutputFile)
		{
			fOutputFile -> Close(); 

			delete fOutputFile; 
		}

	// remember current directory 

	TDirectory * dir = gDirectory; 

	// create a new file 

	fOutputFile = new TFile(fFilename, "RECREATE"); 	

	// change back to the old directory 

	gDirectory = dir; 

	// initialize tree 

	this -> InitializeTree(); 

}

// --------------------------------------------------------- 
	
void BCModelOutput::Fill()
{

	// get output values from model 

	fNParameters = fModel -> GetNParameters(); 
	fProbability_apriori   = fModel -> GetModelAPrioriProbability(); 
	fProbability_aposteriori = fModel -> GetModelAPosterioriProbability(); 
	
	// loop over parameters 

	int nparameters = fModel -> GetNParameters(); 

	for (int i = 0; i < nparameters; ++i)
		{
			BCParameter * parameter = fModel -> GetParameter(i); 

			if (fModel -> GetBestFitParameters().size() > 0)
				fMode_global[i] = fModel -> GetBestFitParameters().at(i); 

			if (fModel -> GetMarginalized(parameter -> GetName()))
				{
					fMode_marginalized[i] = fModel -> GetMarginalized(parameter -> GetName()) -> GetMode(); 
					fMean_marginalized[i] = fModel -> GetMarginalized(parameter -> GetName()) -> GetMean(); 
 					fMedian_marginalized[i] = fModel -> GetMarginalized(parameter -> GetName()) -> GetMedian(); 
 					fQuantile_05[i] = fModel -> GetMarginalized(parameter -> GetName()) -> GetQuantile(0.05); 
 					fQuantile_10[i] = fModel -> GetMarginalized(parameter -> GetName()) -> GetQuantile(0.10); 
 					fQuantile_16[i] = fModel -> GetMarginalized(parameter -> GetName()) -> GetQuantile(0.16); 
 					fQuantile_84[i] = fModel -> GetMarginalized(parameter -> GetName()) -> GetQuantile(0.84); 
 					fQuantile_90[i] = fModel -> GetMarginalized(parameter -> GetName()) -> GetQuantile(0.90); 
 					fQuantile_95[i] = fModel -> GetMarginalized(parameter -> GetName()) -> GetQuantile(0.95); 
				}
		}

	// fill tree 

	fOutputTree -> Fill(); 

	// increase index 

	fIndex++; 

}

// --------------------------------------------------------- 
	
void BCModelOutput::Close()
{

	// remember current directory 

	TDirectory * dir = gDirectory;

	// change to file 

	fOutputFile -> cd(); 

	// write to file 

	fOutputTree -> Write(); 

	// close file 

	fOutputFile -> Close(); 

	// return to old directory 

	gDirectory = dir; 

}

// --------------------------------------------------------- 
	
void BCModelOutput::InitializeTree() 
{

	// create new tree 

	if (fOutputTree)
		delete fOutputTree; 

	fOutputTree = new TTree("BATtree", "BATtree"); 

	// set branch addresses 

	fOutputTree -> Branch("fIndex",                   &fIndex,                   "index/I"); 
	fOutputTree -> Branch("fNParameters",             &fNParameters,             "parameters/I"); 
	fOutputTree -> Branch("fProbability_apriori" ,    &fProbability_apriori,     "apriori probability/D"); 
	fOutputTree -> Branch("fProbability_aposteriori", &fProbability_aposteriori, "aposteriori probability/D"); 
	fOutputTree -> Branch("fMode_global",              fMode_global,             "mode (global) [parameters]/D"); 
	fOutputTree -> Branch("fMode_marginalized",        fMode_marginalized,       "mode (marginalized) [parameters]/D"); 
	fOutputTree -> Branch("fMean_marginalized",        fMean_marginalized,       "mean (marginalized)[parameters]/D"); 
	fOutputTree -> Branch("fMedian_marginalized",      fMedian_marginalized,     "median (marginalized)[parameters]/D"); 
	fOutputTree -> Branch("fQuantile_05" ,             fQuantile_05,             "quantile 5% [parameters]/D"); 
	fOutputTree -> Branch("fQuantile_10" ,             fQuantile_10,             "quantile 10% [parameters]/D"); 
	fOutputTree -> Branch("fQuantile_16" ,             fQuantile_16,             "quantile 16% [parameters]/D"); 
	fOutputTree -> Branch("fQuantile_84" ,             fQuantile_84,             "quantile 84% [parameters]/D"); 
	fOutputTree -> Branch("fQuantile_90" ,             fQuantile_90,             "quantile 90% [parameters]/D"); 
	fOutputTree -> Branch("fQuantile_95" ,             fQuantile_95,             "quantile 95% [parameters]/D"); 

}

// --------------------------------------------------------- 
