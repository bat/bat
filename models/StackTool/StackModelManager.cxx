#include <iostream>
#include <fstream> 

#include "StackModelManager.h"

#include <TROOT.h>
#include <THStack.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TMath.h>

// ---------------------------------------------------------
StackModelManager::StackModelManager()
{
}

// ---------------------------------------------------------
StackModelManager::~StackModelManager()
{}

// ---------------------------------------------------------
int StackModelManager::SetDataHistogram(TH1D hist)
{
	int n = int(fStackModelContainer.size());
	int errcode = 1;

	// loop over stack models
	for (int i = 0; i < n; ++i)
	{
		StackModel * sm = fStackModelContainer.at(i);
		errcode *= sm->SetDataHistogram(hist);
	}

	// return error code
	return errcode;
}

// ---------------------------------------------------------
void StackModelManager::SetFlagFixNorm(bool flag)
{
	int n = int(fStackModelContainer.size());

	// loop over stack models
	for (int i = 0; i < n; ++i)
	{
		StackModel * sm = fStackModelContainer.at(i);
		sm->SetFlagFixNorm(flag);
	}
}

// ---------------------------------------------------------
void StackModelManager::AddStackModel(StackModel * sm)
{
	fStackModelContainer.push_back(sm);

	fNParameters.push_back(sm->GetNParameters());
}

// ---------------------------------------------------------
int StackModelManager::PerformAnalysis()
{
	int n = int(fStackModelContainer.size());
	int errcode = 1;

	// reset
	fChi2.clear();
	fChi2Prob.clear();
	fKSProb.clear();
	fPValue.clear();
	fMaxLike.clear();
	fPosterior.clear();
	fKFactor.clear();
	fNDF.clear();
	fNParameters.clear();

	double normsum = 0.;

	// loop over stack models
	for (int i = 0; i < n; ++i)
	{
		StackModel * sm = fStackModelContainer.at(i);

		// set options
		sm->MCMCSetFlagOrderParameters(false);
		sm->SetOptimizationMethod(BCIntegrate::kOptMinuit);

		// fill basic parameters
		fNDF.push_back(sm->GetNDF());
		fNParameters.push_back(sm->GetNParameters());

		// calculate normalization
		sm->Normalize();
		double norm = sm->GetNormalization();
		fNorm.push_back(norm);
		normsum += norm;

		// run Markov Chains
		sm->MarginalizeAll();

		// find global mode
		sm->FindMode();
		double maxlike = sm->CalculateMaxLike();
		fMaxLike.push_back(maxlike);

		// calculate chi2
		double chi2 = sm->CalculateChi2();
		fChi2.push_back(chi2);

		// calculate chi2prob
		double chi2prob = sm->CalculateChi2Prob();
		fChi2Prob.push_back(chi2prob);

		// do Kolmogorov-Smirnov test
		double probks = sm->CalculateKSProb();
		fKSProb.push_back(probks);

		// calculate p-value
		double pval = sm->CalculatePValue();
		fPValue.push_back(pval);
	}

	// loop over stack models
	for (int i = 0; i < n; ++i)
	{
		// calculate posterior
		fPosterior.push_back(fNorm.at(i) / normsum);

		// calculate k-factor
		std::vector <double> vec;
		for (int j = 0; j < n; ++j)
			vec.push_back(fNorm.at(i) / fNorm.at(j));
		fKFactor.push_back(vec);
	}

	// return error code
	return errcode;
}

// ---------------------------------------------------------
void StackModelManager::PrintResults(const char * filename)
{
	// open file
	std::ofstream ofi(filename);

	// check if file is open
	if(!ofi.is_open())
		{
			std::cerr << "Couldn't open file " << filename <<std::endl;
			return;
		}

	int n = int(fStackModelContainer.size());

	ofi << std::endl;
	ofi << "Goodness-of-fit results " << std::endl;
	ofi << "----------------------- " << std::endl;
	ofi << std::endl;

	ofi << "chi2 / ndf (chiprob): " << std::endl;
	ofi << std::endl;

	// loop over stack models
	for (int i = 0; i < n; ++i)
		ofi << " model " << i << " : " << fChi2.at(i) << " / " << fNDF.at(i) << " (" << fChi2Prob.at(i) << ")" << std::endl;
	ofi << std::endl;

	ofi << "Kolmogorov-Smirnov probability: " << std::endl;
	ofi << std::endl;

	// loop over stack models
	for (int i = 0; i < n; ++i)
		ofi << " model " << i << " : " << fKSProb.at(i) << std::endl;
	ofi << std::endl;

	ofi << "p-value: " << std::endl;
	ofi << std::endl;

	// loop over stack models
	for (int i = 0; i < n; ++i)
		ofi << " model " << i << " : " << fPValue.at(i) << std::endl;

	ofi << std::endl;
	ofi << "Model comparison results " << std::endl;
	ofi << "------------------------ " << std::endl;
	ofi << std::endl;

	ofi << "Maximum Likelihood: " << std::endl;
	ofi << std::endl;
	// loop over stack models
	for (int i = 0; i < n; ++i)
		ofi << " model " << i << " : " << fMaxLike.at(i) << std::endl;
	ofi << std::endl;


	ofi << "posterior probability: " << std::endl;
	ofi << std::endl;

	// loop over stack models
	for (int i = 0; i < n; ++i)
		ofi << " model " << i << " : " << fPosterior.at(i) << std::endl;
	ofi << std::endl;

	ofi << "k-factors: " << std::endl;
	ofi << std::endl;

	// loop over stack models
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < n; ++j)
		{
			if (i == j)
				continue;
			ofi << " model " << i << " vs " << j << " : " << (fKFactor.at(i)).at(j) << std::endl;
		}
	ofi << std::endl;

}

// ---------------------------------------------------------
