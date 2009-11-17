#include <iostream>

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
void StackModelManager::PrintResults()
{
	int n = int(fStackModelContainer.size());

	std::cout << std::endl;
	std::cout << "Goodness-of-fit results " << std::endl;
	std::cout << "----------------------- " << std::endl;
	std::cout << std::endl;

	std::cout << "chi2 / ndf (chiprob): " << std::endl;
	std::cout << std::endl;

	// loop over stack models
	for (int i = 0; i < n; ++i)
		std::cout << " model " << i << " : " << fChi2.at(i) << " / " << fNDF.at(i) << " (" << fChi2Prob.at(i) << ")" << std::endl;
	std::cout << std::endl;

	std::cout << "Kolmogorov-Smirnov probability: " << std::endl;
	std::cout << std::endl;

	// loop over stack models
	for (int i = 0; i < n; ++i)
		std::cout << " model " << i << " : " << fKSProb.at(i) << std::endl;
	std::cout << std::endl;

	std::cout << "p-value: " << std::endl;
	std::cout << std::endl;

	// loop over stack models
	for (int i = 0; i < n; ++i)
		std::cout << " model " << i << " : " << fPValue.at(i) << std::endl;

	std::cout << std::endl;
	std::cout << "Model comparison results " << std::endl;
	std::cout << "------------------------ " << std::endl;
	std::cout << std::endl;

	std::cout << "Maximum Likelihood: " << std::endl;
	std::cout << std::endl;
	// loop over stack models
	for (int i = 0; i < n; ++i)
		std::cout << " model " << i << " : " << fMaxLike.at(i) << std::endl;
	std::cout << std::endl;


	std::cout << "posterior probability: " << std::endl;
	std::cout << std::endl;

	// loop over stack models
	for (int i = 0; i < n; ++i)
		std::cout << " model " << i << " : " << fPosterior.at(i) << std::endl;
	std::cout << std::endl;

	std::cout << "k-factors: " << std::endl;
	std::cout << std::endl;

	// loop over stack models
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < n; ++j)
		{
			if (i == j)
				continue;
			std::cout << " model " << i << " vs " << j << " : " << (fKFactor.at(i)).at(j) << std::endl;
		}
	std::cout << std::endl;

}

// ---------------------------------------------------------
