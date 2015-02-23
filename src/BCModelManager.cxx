/*
 * Copyright (C) 2007-2013, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

// ---------------------------------------------------------

#include "BCModelManager.h"

#include "BCDataPoint.h"
#include "BCLog.h"

#include <TString.h>

#include <fstream>
#include <iostream>

// ---------------------------------------------------------
BCModelManager::BCModelManager()
	: fDataSet(0)
{
}

// ---------------------------------------------------------

BCModelManager::~BCModelManager() {
	delete fDataSet;
}

// ---------------------------------------------------------

BCModelManager::BCModelManager(const BCModelManager & other)
	: fModels(other.fModels)
	, fDataSet(other.fDataSet)
{
}

// ---------------------------------------------------------

BCModelManager & BCModelManager::operator = (const BCModelManager & rhs) {
	// copy pointers only
	fModels  = rhs.fModels;
	fDataSet = rhs.fDataSet;
	return *this;
}

// ---------------------------------------------------------

void BCModelManager::SetDataSet(BCDataSet * dataset) {
	// set data set
	fDataSet = dataset;

	// set data set of all models in the manager
	for (unsigned int i = 0; i < GetNModels(); ++i)
		GetModel(i) -> SetDataSet(fDataSet);
}

// ---------------------------------------------------------
void BCModelManager::AddModel(BCModel * model, double probability) {
	// set data set
	model -> SetDataSet(fDataSet);

	// fill model into container
	fModels.push_back(model);
	// set probabilities
	fAPrioriProbability.push_back(probability);
	fAPosterioriProbability.push_back(-1);
}

// ---------------------------------------------------------
void BCModelManager::MCMCSetPrecision(BCEngineMCMC::Precision precision) {
	for (unsigned i = 0; i < GetNModels(); ++i)
		GetModel(i) -> MCMCSetPrecision(precision);
}

// ---------------------------------------------------------
void BCModelManager::SetNIterationsMax(int niterations) {
	for (unsigned i = 0; i < GetNModels(); ++i)
		GetModel(i) -> SetNIterationsMax(niterations);
}

// ---------------------------------------------------------
void BCModelManager::SetNIterationsMin(int niterations) {
	for (unsigned i = 0; i < GetNModels(); ++i)
		GetModel(i) -> SetNIterationsMin(niterations);
}

// ---------------------------------------------------------
void BCModelManager::SetNIterationsPrecisionCheck(int niterations) {
	for (unsigned i = 0; i < GetNModels(); ++i)
		GetModel(i) -> SetNIterationsPrecisionCheck(niterations);
}

// ---------------------------------------------------------
void BCModelManager::SetIntegrationMethod(BCIntegrate::BCIntegrationMethod method) {
	for (unsigned i = 0; i < GetNModels(); ++i)
		GetModel(i) -> SetIntegrationMethod(method);
}

// ---------------------------------------------------------
void BCModelManager::SetMarginalizationMethod(BCIntegrate::BCMarginalizationMethod method) {
	for (unsigned i = 0; i < GetNModels(); ++i)
		GetModel(i) -> SetMarginalizationMethod(method);
}

// ---------------------------------------------------------
void BCModelManager::SetOptimizationMethod(BCIntegrate::BCOptimizationMethod method) {
	for (unsigned i = 0; i < GetNModels(); ++i)
		GetModel(i) -> SetOptimizationMethod(method);
}

// ---------------------------------------------------------
void BCModelManager::SetRelativePrecision(double relprecision) {
	for (unsigned i = 0; i < GetNModels(); ++i)
		GetModel(i) -> SetRelativePrecision(relprecision);
}

// ---------------------------------------------------------
void BCModelManager::SetAbsolutePrecision(double absprecision) {
	for (unsigned i = 0; i < GetNModels(); ++i)
		GetModel(i) -> SetAbsolutePrecision(absprecision);
}

// ---------------------------------------------------------
void BCModelManager::SetNbins(unsigned int n) {
	for (unsigned i = 0; i < GetNModels(); ++i)
		GetModel(i) -> SetNbins(n);
}

// ---------------------------------------------------------
void BCModelManager::SetNChains(unsigned int n) {
	for (unsigned i = 0; i < GetNModels(); ++i)
		GetModel(i) -> MCMCSetNChains(n);
}

// ---------------------------------------------------------
void BCModelManager::Integrate() {
	// initialize likelihood norm
	double normalization = 0.0;

	BCLog::OutSummary("Running normalization of all models.");

	for (unsigned int i = 0; i < GetNModels(); ++i) {
		fModels[i] -> Integrate();

		// add to total normalization
		normalization += fModels[i]->GetIntegral() * fAPrioriProbability[i];
	}

	// set model a posteriori probabilities
	for (unsigned int i = 0; i < GetNModels(); ++i)
		fAPosterioriProbability[i] = (fModels[i]->GetIntegral() * fAPrioriProbability[i]) / normalization;
}

// ---------------------------------------------------------
double BCModelManager::BayesFactor(unsigned imodel1, unsigned imodel2) {
	if (imodel1 >= fModels.size() or imodel2 >= fModels.size())
		return -1;

	// Bayes Factors are the likelihoods integrated over the parameters
	// Is this equal to the posteriors?
	//    NOOOO
	// But it is equal to normalization factors.
	
	// check model 1
	const double norm1 = fModels[imodel1] -> GetIntegral();
	if (norm1<0.) {
		BCLog::OutError(Form("BCModelManager::BayesFactor : Model %s (index %d) not normalized. Cannot calculate Bayes factor.",fModels[imodel1]->GetName().data(),imodel1));
		return -1.;
	}

	// check model 2
	const double norm2 = fModels[imodel2] -> GetIntegral();
	if (norm2<0) {
		BCLog::OutError(Form("BCModelManager::BayesFactor : Model %s (index %d) not normalized. Cannot calculate Bayes factor.",fModels[imodel2]->GetName().data(),imodel2));
		return -1.;
	}

	// denominator cannot be zero
	if (norm2 < std::numeric_limits<double>::epsilon() and norm1 >= std::numeric_limits<double>::epsilon()) {
		BCLog::OutError(Form("BCModelManager::BayesFactor : Model %s (index %d) has ZERO probability. Bayes factor is infinite.",fModels[imodel2]->GetName().data(),imodel2));
		return -1.;
	}

	// denominator cannot be zero unless also numerator is zero
	if(norm2 < std::numeric_limits<double>::epsilon() and norm1 < std::numeric_limits<double>::epsilon()) {
		BCLog::OutWarning(Form("BCModelManager::BayesFactor : Models %s and %s have ZERO probability. Bayes factor is unknown. Returning 1.",fModels[imodel2]->GetName().data(),fModels[imodel1]->GetName().data()));
		return 1.;
	}
	
	// now calculate the factor
	return norm1/norm2;
}

// ---------------------------------------------------------
void BCModelManager::FindMode() {
	for (unsigned i = 0; i < GetNModels(); ++i)
		GetModel(i) -> FindMode();
}

// ---------------------------------------------------------
void BCModelManager::MarginalizeAll() {
	// marginalizes all models registered
	for (unsigned i = 0; i < GetNModels(); ++i)
		GetModel(i) -> MarginalizeAll();
}

// ---------------------------------------------------------
void BCModelManager::WriteMarkovChain(bool flag) {
	for (unsigned i = 0; i < GetNModels(); ++i)
		GetModel(i) -> WriteMarkovChain(flag);
}

// ---------------------------------------------------------
void BCModelManager::WriteMarkovChain(std::string prefix, std::string option) {
	for (unsigned i = 0; i < GetNModels(); ++i)
		GetModel(i) -> WriteMarkovChain(prefix+GetModel(i)->GetSafeName()+".root",option);
}

// ---------------------------------------------------------
// void BCModelManager::CalculatePValue(bool flag_histogram)
// {
// 	for (unsigned i = 0; i < GetNModels(); ++i)
// 		GetModel(i) -> CalculatePValue(GetModel(i)->GetBestFitParameters(), flag_histogram);
// }

// ---------------------------------------------------------
void BCModelManager::PrintSummary(const char * file) {
  std::ofstream out;
  std::streambuf * old_buffer = 0;
	
	if(file) {
		out.open(file);
		if (!out.is_open()) {
			std::cerr<<"Couldn't open file "<<file<<std::endl;
			return;
		}
		old_buffer = std::cout.rdbuf(out.rdbuf());
	}

	// model summary
	std::cout << std::endl
					  << "======================================" << std::endl
					  << " Summary" << std::endl
					  << "======================================" << std::endl
					  << std::endl
					  << " Number of models               : " << fModels.size() << std::endl
					  << std::endl
					  << " - Models:" << std::endl;

	for (unsigned i = 0; i < fModels.size(); ++i)
		fModels[i] -> PrintSummary();
	
	// data summary
	std::cout << " - Data:" << std::endl
					  << std::endl
					  << "     Number of entries: " << fDataSet->GetNDataPoints() << std::endl << std::endl;

	std::cout << "======================================" << std::endl
					  << " Model comparison" << std::endl << std::endl;

	// probability summary
	std::cout <<" - A priori probabilities:" << std::endl << std::endl;
	
	for (unsigned i=0; i<fModels.size(); ++i)
		std::cout << "     p(" << fModels[i]->GetName()	<< ") = " << fAPrioriProbability[i] << std::endl;

	std::cout << std::endl << " - A posteriori probabilities:" << std::endl << std::endl;

	for (unsigned i = 0; i < fModels.size(); ++i)
		std::cout << "     p(" << fModels[i]->GetName()	<< " | data) = " << fAPosterioriProbability[i] << std::endl;

	std::cout << std::endl << "======================================" << std::endl << std::endl;
	
	if (file)
		std::cout.rdbuf(old_buffer);
}

// ---------------------------------------------------------

void BCModelManager::PrintModelComparisonSummary(const char * file)
{
	std::ofstream out;
	std::streambuf * old_buffer = 0;
	
	if(file) {
		out.open(file);
		if (!out.is_open()) {
			std::cerr<<"Couldn't open file "<<file<<std::endl;
			return;
		}
		old_buffer = std::cout.rdbuf(out.rdbuf());
	}

	// model summary
	std::cout << std::endl
					  << "===========================================" << std::endl
					  << " Model Comparison Summary" << std::endl
					  << "===========================================" << std::endl
					  << std::endl
					  << " Number of models               : "<< fModels.size() << std::endl << std::endl;
	
	// probability summary
	std::cout << " - A priori probabilities:" << std::endl << std::endl;

	for (unsigned i=0; i<fModels.size(); ++i)
		std::cout << "     p(" << fModels[i]->GetName() << ") = " << fAPrioriProbability[i] << std::endl;

	std::cout << std::endl << " - A posteriori probabilities:" << std::endl << std::endl;

	for (unsigned i = 0; i < fModels.size(); ++i)
		std::cout << "     p(" << fModels[i]->GetName() << " | data) = " << fAPosterioriProbability[i] << std::endl;

	// Bayes factors summary
	std::cout << std::endl << " - Bayes factors:" << std::endl << std::endl;
	for (unsigned i = 0; i < fModels.size(); ++i)
		for (unsigned j = i+1; j < fModels.size(); ++j)
			std::cout<<"     K = p(data | " << fModels[i]->GetName() << ") / p(data | " << fModels[j]->GetName() << ") = " << BayesFactor(i,j) << std::endl;
	
	// // p-values summary
	// std::cout	<< std::endl << " - p-values:" << std::endl << std::endl;

	// for (unsigned i = 0; i < fModels.size(); ++i) {
	// 	double p = fModels[i] -> GetPValue();
	// 	std::cout << "     " << fModels[i]->GetName();
	// 	if(p>=0.)
	// 		std::cout << ":  p-value = " << p << std::endl;
	// 	else
	// 		std::cout << ":  p-value not calculated" << std::endl;
	// }

	std::cout << std::endl << "===========================================" << std::endl << std::endl;

	if (file)
		std::cout.rdbuf(old_buffer);
}

// ---------------------------------------------------------

void BCModelManager::PrintResults()
{
	// print summary of all models
	for (unsigned i = 0; i < GetNModels(); ++i)
		GetModel(i) -> PrintResults(Form("%s.txt", GetModel(i)->GetSafeName().data()));
}
