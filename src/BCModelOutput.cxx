#include "BCModelOutput.h" 

#include "TDirectory.h"

// --------------------------------------------------------- 

BCModelOutput::BCModelOutput() 
{

	fIndex = 0; 

	fOutputFile = 0; 
	fMarkovChainTree = 0; 
	fAnalysisTree = 0; 
	fModel = 0; 
	fOutputFile = 0; 

}

// --------------------------------------------------------- 

BCModelOutput::BCModelOutput(BCModel * model, const char * filename) 
{

	// remeber current directory 

	TDirectory * dir = gDirectory; 

	fFilename = filename; 

	// create a new file 

	fOutputFile = new TFile(fFilename, "RECREATE"); 

	BCModelOutput(); 

	fModel = model;

	// initialize trees 

	this -> InitializeAnalysisTree(); 

	this -> InitializeMarkovChainTree(); 

	// change again to the old directory 

	gDirectory = dir; 

}

// --------------------------------------------------------- 

BCModelOutput::~BCModelOutput() 
{

	if (fOutputFile)
		{
			fOutputFile -> Close(); 

			delete fOutputFile; 
		} 

	if (fAnalysisTree)
		delete fAnalysisTree; 

	if (fMarkovChainTree)
		delete fMarkovChainTree; 

	if (fOutputFile)
		delete fOutputFile; 

	fModel = 0; 

}

// --------------------------------------------------------- 

BCModelOutput::BCModelOutput(const BCModelOutput & modeloutput)
{

	modeloutput.Copy(*this); 

}

// --------------------------------------------------------- 

BCModelOutput & BCModelOutput::operator = (const BCModelOutput & modeloutput)
{

	if (this != &modeloutput) 
		modeloutput.Copy(* this); 

	return * this; 

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

	this -> InitializeAnalysisTree(); 

}

// --------------------------------------------------------- 

void BCModelOutput::WriteMarkovChain(bool flag)
{

	if (fModel)
		fModel -> WriteMarkovChain(flag); 

}

// --------------------------------------------------------- 
	
void BCModelOutput::FillAnalysisTree()
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

			if (fModel -> GetMarginalized(parameter -> GetName().data()))
				{
					fMode_marginalized[i] = fModel -> GetMarginalized(parameter -> GetName().data()) -> GetMode(); 
					fMean_marginalized[i] = fModel -> GetMarginalized(parameter -> GetName().data()) -> GetMean(); 
 					fMedian_marginalized[i] = fModel -> GetMarginalized(parameter -> GetName().data()) -> GetMedian(); 
 					fQuantile_05[i] = fModel -> GetMarginalized(parameter -> GetName().data()) -> GetQuantile(0.05); 
 					fQuantile_10[i] = fModel -> GetMarginalized(parameter -> GetName().data()) -> GetQuantile(0.10); 
 					fQuantile_16[i] = fModel -> GetMarginalized(parameter -> GetName().data()) -> GetQuantile(0.16); 
 					fQuantile_84[i] = fModel -> GetMarginalized(parameter -> GetName().data()) -> GetQuantile(0.84); 
 					fQuantile_90[i] = fModel -> GetMarginalized(parameter -> GetName().data()) -> GetQuantile(0.90); 
 					fQuantile_95[i] = fModel -> GetMarginalized(parameter -> GetName().data()) -> GetQuantile(0.95); 
				}
		}

	// fill tree 

	fAnalysisTree -> Fill(); 

	// increase index 

	fIndex++; 

}

// --------------------------------------------------------- 

void BCModelOutput::WriteMarginalizedDistributions() 
{

	// remember current directory 

	TDirectory * dir = gDirectory;

	// change to file 

	fOutputFile -> cd(); 

	int nparameters = fModel -> GetNParameters(); 

	for (int i = 0; i < nparameters; ++i)
		fModel -> GetMarginalized(fModel -> GetParameter(i)) -> GetHistogram() -> Write(); 

	if (nparameters > 1) 
		{
			for (int i = 0; i < nparameters - 1; ++i)
				for (int j = i + 1; j < nparameters; ++j) 
					fModel -> GetMarginalized(fModel -> GetParameter(i), 
																		fModel -> GetParameter(j)) -> GetHistogram() -> Write(); 
		}

	// return to old directory 

	gDirectory = dir; 

}

// --------------------------------------------------------- 
	
void BCModelOutput::Close()
{

	// remember current directory 

	TDirectory * dir = gDirectory;

	// change to file 

	fOutputFile -> cd(); 

	// write analysis tree to file 

	if (fAnalysisTree -> GetEntries() > 0) 
		fAnalysisTree -> Write(); 

	// write markov chain tree to file 

	if (fMarkovChainTree -> GetEntries() > 0) 
		fMarkovChainTree -> Write(); 

	// close file 

	fOutputFile -> Close(); 

	// return to old directory 

	gDirectory = dir; 

}

// --------------------------------------------------------- 
	
void BCModelOutput::InitializeAnalysisTree() 
{

	// create new tree 

	if (fAnalysisTree)
		delete fAnalysisTree; 

	fAnalysisTree = new TTree("AnalysisTree", "AnalysisTree"); 

	// set branch addresses 

	fAnalysisTree -> Branch("fIndex",                   &fIndex,                   "index/I"); 
	fAnalysisTree -> Branch("fNParameters",             &fNParameters,             "parameters/I"); 
	fAnalysisTree -> Branch("fProbability_apriori" ,    &fProbability_apriori,     "apriori probability/D"); 
	fAnalysisTree -> Branch("fProbability_aposteriori", &fProbability_aposteriori, "aposteriori probability/D"); 
	fAnalysisTree -> Branch("fMode_global",              fMode_global,             "mode (global) [parameters]/D"); 
	fAnalysisTree -> Branch("fMode_marginalized",        fMode_marginalized,       "mode (marginalized) [parameters]/D"); 
	fAnalysisTree -> Branch("fMean_marginalized",        fMean_marginalized,       "mean (marginalized)[parameters]/D"); 
	fAnalysisTree -> Branch("fMedian_marginalized",      fMedian_marginalized,     "median (marginalized)[parameters]/D"); 
	fAnalysisTree -> Branch("fQuantile_05" ,             fQuantile_05,             "quantile 5% [parameters]/D"); 
	fAnalysisTree -> Branch("fQuantile_10" ,             fQuantile_10,             "quantile 10% [parameters]/D"); 
	fAnalysisTree -> Branch("fQuantile_16" ,             fQuantile_16,             "quantile 16% [parameters]/D"); 
	fAnalysisTree -> Branch("fQuantile_84" ,             fQuantile_84,             "quantile 84% [parameters]/D"); 
	fAnalysisTree -> Branch("fQuantile_90" ,             fQuantile_90,             "quantile 90% [parameters]/D"); 
	fAnalysisTree -> Branch("fQuantile_95" ,             fQuantile_95,             "quantile 95% [parameters]/D"); 

}


// --------------------------------------------------------- 
	
void BCModelOutput::InitializeMarkovChainTree() 
{

	// create new tree 

	if (fMarkovChainTree)
		delete fMarkovChainTree; 

	fMarkovChainTree = new TTree("MarkovChainTree", "MarkovChainTree"); 

	// connect pointer to parameter vectors 

	fParameters = fModel -> GetMarkovChainPoint(); 

	fIteration = fModel -> GetMCMCIteration(); 

	fLogLikelihood = fModel -> GetMarkovChainValue(); 

	fMarkovChainTree -> Branch("fIteration",      fIteration,   "iteration/I"); 
	fMarkovChainTree -> Branch("fNParameters",   &fNParameters, "parameters/I"); 
	fMarkovChainTree -> Branch("fLogLikelihood",  fLogLikelihood, "log (likelihood)/D"); 

	// loop over all parameters 

	fNParameters = fModel -> GetNParameters(); 

	for (int i = 0; i < fNParameters; ++i)
		{
			fMarkovChainTree -> Branch(Form("fParameter%i", i), 
																 &(*fParameters)[i], 
																 Form("parameter %i/D", i)); 
		}

	fModel -> SetMarkovChainTree(fMarkovChainTree); 

}

// --------------------------------------------------------- 

void BCModelOutput::Copy(BCModelOutput & modeloutput) const 
{

	// don't copy the content 

	modeloutput.fModel           = this -> fModel; 
	modeloutput.fAnalysisTree    = this -> fAnalysisTree; 
	modeloutput.fMarkovChainTree = this -> fMarkovChainTree; 

}

// --------------------------------------------------------- 
