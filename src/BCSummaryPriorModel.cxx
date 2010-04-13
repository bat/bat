#include <TH1D.h>

#include <BAT/BCLog.h>
#include <BAT/BCH1D.h>

#include <BAT/BCSummaryPriorModel.h>

// ---------------------------------------------------------
BCSummaryPriorModel::BCSummaryPriorModel()
	: BCModel()
	, fTestModel(0)
{
}

// ---------------------------------------------------------
BCSummaryPriorModel::BCSummaryPriorModel(const char * name)
	: BCModel(name)
	, fTestModel(0)
{
}

// ---------------------------------------------------------
BCSummaryPriorModel::~BCSummaryPriorModel()
{}

// ---------------------------------------------------------
double BCSummaryPriorModel::LogLikelihood(std::vector <double> parameters)
{
	return fTestModel->LogAPrioriProbability(parameters);
}

// ---------------------------------------------------------
double BCSummaryPriorModel::LogAPrioriProbability(std::vector <double> parameters)
{
	return 0;
}

// ---------------------------------------------------------
int BCSummaryPriorModel::PerformAnalysis()
{
	// check if model is set
	if (!fTestModel) {
		BCLog::OutError("BCSummaryPriorModel::PerformAnalysis : Model not defined.");
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

