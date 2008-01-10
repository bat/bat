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

	this -> InitializeMarkovChainTrees(); 

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

	for (int i = 0; i < fModel -> MCMCGetNChains(); ++i)
		if (fMarkovChainTrees[i] -> GetEntries() > 0) 
			fMarkovChainTrees[i] -> Write(); 

	// write control plots to file 

	if (fModel -> MCMCGetH1RValue())
	  fModel -> MCMCGetH1RValue() -> Write(); 

	if (fModel -> MCMCGetH1Efficiency())
	  fModel -> MCMCGetH1Efficiency() -> Write(); 

	// close file 

	fOutputFile -> Close(); 

	// return to old directory 

	gDirectory = dir; 

}

// --------------------------------------------------------- 
	
void BCModelOutput::InitializeAnalysisTree() 
{

	// create new tree 

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
	
void BCModelOutput::InitializeMarkovChainTrees() 
{

	// create new tree 

	fMarkovChainTrees.clear(); 

	for (int i = 0; i < fModel -> MCMCGetNChains(); ++i)
		{
			TTree * tree = new TTree(Form("MarkovChainTree_%i", i), "MarkovChainTree"); 
			fMarkovChainTrees.push_back(tree); 
		}

// 	// connect pointer to parameter vectors 

	fParameters    = fModel -> MCMCGetP2x(); 
	fIteration     = fModel -> MCMCGetP2NIterations(); 
	fLogLikelihood = fModel -> MCMCGetP2LogProbx(); 

	fNParameters = fModel -> MCMCGetNParameters(); 

	for (int i = 0; i < fModel -> MCMCGetNChains(); ++i)
		{
			fMarkovChainTrees[i] -> Branch("fIteration",      &(*fIteration)[i],   "iteration/I"); 
			fMarkovChainTrees[i] -> Branch("fNParameters",    &fNParameters,       "parameters/I"); 
			fMarkovChainTrees[i] -> Branch("fLogLikelihood",  &(*fLogLikelihood)[i], "log (likelihood)/D"); 
		}

	// loop over all parameters 
	
	for (int i = 0; i < fModel -> MCMCGetNChains(); ++i)
		for (int j = 0; j <  fModel -> MCMCGetNParameters(); ++j)
		{
			fMarkovChainTrees[i] -> Branch(Form("fParameter%i", j), 
																		 &(*fParameters)[i * fModel -> MCMCGetNParameters() + j], 
																		 Form("parameter %i/D", j)); 
		}

	fModel -> MCMCSetMarkovChainTrees(fMarkovChainTrees); 

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
