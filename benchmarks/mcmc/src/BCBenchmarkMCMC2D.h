/*!
 * \class BCBenchmarkMCMC2D
 * \brief A class to check 2D MCMC properties
 * \author Daniel Kollar
 * \author Kevin Kr&ouml;ninger
 * \author Jing Liu
 * \version 1.0
 * \date 01.09.2009
 * \detail This class create MCMCs according to a known 2D distribution 
 * and check the statistical properties of the created chains, e.g.
 * x2, moments, quantiles and fluctuations, etc. as functions of 
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

#ifndef __BCBENCHMARKMCMC2D__H
#define __BCBENCHMARKMCMC2D__H

#include <TF2.h>
#include <TH2F.h>

#include "BAT/BCModel.h"
#include "BAT/BCModelOutput.h"

//=============================================================================

class BCBenchmarkMCMC2D : public BCModel, public BCModelOutput
{
	public:

		BCBenchmarkMCMC2D(
				TF2* testFunction = NULL,
				const char* outputFile = "MCMCtest2d.root",
				const char* modelName = "BenchmarkMCMC2D");
		~BCBenchmarkMCMC2D();

		// inherited methods
		double LogAPrioriProbability(std::vector <double> parameters)
		{return 0;}

		double LogLikelihood(std::vector <double> parameters)
		{return log(fTestFunction->Eval(parameters[0],parameters[1]));}

		// own methods
		void ProcessMCTrees();

		void PerformLagsTest();
		void PerformIterationsTest();

		void WriteResults();

	private:

		void ProcessMCTree(int chainID=0);
		void Chi2vsLagsOfChain(int chainID=0);
		void Chi2vsIterOfChain(int chainID=0);
		
		void KolmogorovVsLagsOfChain(int chainID=0);
		void KolmogorovVsIterOfChain(int chainID=0);

		int fNbinx, fNbiny;
		double fXmin, fXmax, fYmin, fYmax;

		static const int fMaxChains = 5;
		static const int fMaxLags = 31; // lag = 1,2,...,30
		static const int fMax10thOfIters = 11; // iterations/10*(1,2,...10)

		TF2* fTestFunction;
		TF2* fFitFunction;

		TH2F* fHistXYLags[fMaxChains][fMaxLags];
		TH2F* fHistXYIter[fMaxChains][fMax10thOfIters];

		TH1F* fHChi2vsLags[fMaxChains];
		TH1F* fHChi2vsIter[fMaxChains];

		TH1F* fHKolmogorovProbVsLags[fMaxChains];
		TH1F* fHKolmogorovProbVsIter[fMaxChains];
};

#endif
