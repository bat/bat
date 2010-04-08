#include <TH1D.h>

#include <BAT/BCLog.h>
#include <BAT/BCH1D.h>

#include <BAT/BCSPriorModel.h>

// ---------------------------------------------------------
BCSPriorModel::BCSPriorModel()
	: BCModel()
	, fTestModel(0)
{
}

// ---------------------------------------------------------
BCSPriorModel::BCSPriorModel(const char * name)
	: BCModel(name)
	, fTestModel(0)
{
}

// ---------------------------------------------------------
BCSPriorModel::~BCSPriorModel()
{}

// ---------------------------------------------------------
double BCSPriorModel::LogLikelihood(std::vector <double> parameters)
{
	return fTestModel->LogAPrioriProbability(parameters);
}

// ---------------------------------------------------------
double BCSPriorModel::LogAPrioriProbability(std::vector <double> parameters)
{
	return 0;
}

// ---------------------------------------------------------
int BCSPriorModel::PerformAnalysis()
{
	// check if model is set
	if (!fTestModel) {
		BCLog::OutError("BCSPriorModel::PerformAnalysis. Model not defined.");
		return 0;
	}

	// copy parameters
	int npar = fTestModel->GetNParameters();
	for (int i = 0; i < npar; ++i) {
		BCParameter* par = fTestModel->GetParameter(i);
		AddParameter(par);
	}

	for (int i = 0; i < npar; ++i) {
		int nbins = fTestModel->GetMarginalized(fTestModel->GetParameter(i))->GetHistogram()->GetNbinsX();
		SetNbins((fTestModel->GetParameter(i)->GetName()).c_str(), nbins);
	}

	MCMCSetNLag(10);
	MCMCSetNIterationsRun( fTestModel->MCMCGetNIterationsRun() );
	MarginalizeAll();
	FindModeMinuit( GetBestFitParameters() );

	// no error
	return 1;
}

// ---------------------------------------------------------

