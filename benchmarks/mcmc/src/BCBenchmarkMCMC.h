/*!
 * \class BCBenchmarkMCMC
 * \brief A class to check MCMC properties
 * \author Daniel Kollar
 * \author Kevin Kr&ouml;ninger
 * \author Jing Liu
 * \version 1.0
 * \date 01.09.2009
 * \detail This class create MCMCs according to a known distribution
 * and check the statistical properties of the created chains, e.g.
 * x2, moments, quantiles and fluctuations, etc as functions of
 * different lags and iterations.
 */

/*
 * Copyright (C) 2009,
 * Daniel Kollar, Kevin Kroeninger and Jing Liu.
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

//=============================================================================

#ifndef __BCBENCHMARKMCMC__H
#define __BCBENCHMARKMCMC__H

#include <TF1.h>
#include <TH1F.h>

#include "BAT/BCModel.h"
#include "BAT/BCModelOutput.h"

//=============================================================================

class BCBenchmarkMCMC : public BCModel, public BCModelOutput
{
 public:
	
	BCBenchmarkMCMC(
									TF1* testFunction = NULL,
									const char* outputFile = "MCMCtest.root",
									const char* modelName = "BenchmarkMCMC"
									);
	~BCBenchmarkMCMC();

	// inherited methods
	double LogAPrioriProbability(std::vector<double> parameters)
	{return 0;}

	double LogLikelihood(std::vector<double> parameters)
	{return log(fTestFunction->Eval(parameters[0]));}

	// own methods
	void ProcessMCTrees();

	void PerformLagsTest();
	void PerformIterationsTest();

	void WriteResults();
	void PrintComparison(const char* filename);

 private:

	void ProcessMCTree(int chainID=0);
	void Chi2vsLagsOfChain(int chainID=0);
	void Chi2vsIterOfChain(int chainID=0);

	void KolmogorovVsLagsOfChain(int chainID=0);
	void KolmogorovVsIterOfChain(int chainID=0);

	int fNbinx;
	double fXmin, fXmax;

	static const int fMaxChains = 5;
	static const int fMaxLags = 31; // lag = 1,2,...,30
	static const int fMax10thOfIters = 11; // iterations/10*(1,2,...10)

	TF1* fTestFunction;
	TF1* fFitFunction;

	TH1F* fHistXLags[fMaxChains][fMaxLags];
	TH1F* fHistXIter[fMaxChains][fMax10thOfIters];

	TH1F* fHChi2vsLags[fMaxChains];
	TH1F* fHChi2vsIter[fMaxChains];

	TH1F* fHKolmogorovProbVsLags[fMaxChains];
	TH1F* fHKolmogorovProbVsIter[fMaxChains];

};

#endif
