#include <BAT/BCLog.h> 
#include <BAT/BCH1D.h> 

#include <TH1D.h>

#include "PriorModel.h"

// ---------------------------------------------------------
PriorModel::PriorModel() : BCModel()
{
	fTestModel = 0; 
};

// ---------------------------------------------------------
PriorModel::PriorModel(const char * name) : BCModel(name)
{
	fTestModel = 0; 
};

// ---------------------------------------------------------
PriorModel::~PriorModel()
{}; 

// ---------------------------------------------------------
double PriorModel::LogLikelihood(std::vector <double> parameters)
{
	double logprob = fTestModel->LogAPrioriProbability(parameters); 

	return logprob;
}

// ---------------------------------------------------------
double PriorModel::LogAPrioriProbability(std::vector <double> parameters)
{
	return 0;
}

// ---------------------------------------------------------
int PriorModel::PerformAnalysis()
{
	// check if model is set
	if (!fTestModel) {
		BCLog::OutError("PriorModel::PerformAnalysis. Model not defined.");
		return 0; 
	}

	// copy parameters
	int npar = fTestModel->GetNParameters(); 
	for (int i = 0; i < npar; ++i) {
		BCParameter* par = fTestModel->GetParameter(i); 
		this->AddParameter(par); 
	}

	for (int i = 0; i < npar; ++i) {
		int nbins = fTestModel->GetMarginalized(fTestModel->GetParameter(i))->GetHistogram()->GetNbinsX(); 
		this->SetNbins((fTestModel->GetParameter(i)->GetName()).c_str(), nbins);
	}

	this->MCMCSetNLag(10); 
	this->MCMCSetNIterationsRun( fTestModel->MCMCGetNIterationsRun() ); 
	this->MarginalizeAll(); 
	this->FindModeMinuit( this->GetBestFitParameters() ); 

	// no error
	return 1;
}

// ---------------------------------------------------------

