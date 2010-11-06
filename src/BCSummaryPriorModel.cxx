/*
 * Copyright (C) 2008-2010, Daniel Kollar and Kevin Kroeninger.
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

// ---------------------------------------------------------

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
void BCSummaryPriorModel::SetModel(BCModel * model)
{
	fTestModel = model;

	// copy parameters
	int npar = fTestModel->GetNParameters();
	for (int i = 0; i < npar; ++i) {
		BCParameter * par = fTestModel->GetParameter(i);
		AddParameter(par);
	}

	// set default histogram binning to the one of the original model
	for (int i = 0; i < npar; ++i) {
		int nbins = fTestModel->GetMarginalized( fTestModel->GetParameter(i) )->GetHistogram()->GetNbinsX();
		SetNbins( (fTestModel->GetParameter(i)->GetName()).c_str(), nbins);
	}

	// set large default lag
	//	MCMCSetNLag(10);

	// set default MCMC setup to the one of the original model
	MCMCSetNChains( fTestModel->MCMCGetNChains() );
	MCMCSetNLag( fTestModel->MCMCGetNLag() );
	MCMCSetNIterationsMax( fTestModel->MCMCGetNIterationsMax() );
	MCMCSetNIterationsRun( fTestModel->MCMCGetNIterationsRun() );
	MCMCSetNIterationsPreRunMin( fTestModel->MCMCGetNIterationsPreRunMin() ); 
	MCMCSetNIterationsUpdate( fTestModel->MCMCGetNIterationsUpdate() );
	MCMCSetNIterationsUpdateMax( fTestModel->MCMCGetNIterationsUpdateMax() );
	MCMCSetRValueCriterion( fTestModel->MCMCGetRValueCriterion() );
	MCMCSetRValueParametersCriterion( fTestModel->MCMCGetRValueParametersCriterion() );

	/*
		fMCMCNChains              = 5;
		fMCMCNLag                 = 1;
		fMCMCNIterationsMax       = 100000;
		fMCMCNIterationsRun       = 100000;
		fMCMCNIterationsPreRunMin = 100;
		fMCMCNIterationsUpdate    = 1000;
		fMCMCNIterationsUpdateMax = 10000;
		fMCMCRValueCriterion      = 0.1;
		fMCMCRValueParametersCriterion = 0.1;
		fMCMCRValue               = 100;
	*/


	//	MCMCSetNChains( fTestModel->MCMCGetNChains() );
	//	MCMCSetNIterationsRun( fTestModel->MCMCGetNIterationsRun() );

}

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

