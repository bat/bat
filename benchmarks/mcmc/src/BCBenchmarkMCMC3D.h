/*!
 * \class BCBenchmarkMCMC3D
 * \brief A class to check 3D MCMC properties
 * \author Daniel Kollar
 * \author Kevin Kr&ouml;ninger
 * \author Jing Liu
 * \author Carsten Brachem
 * \version 1.0
 * \date 16.09.2009
 * \detail This class create MCMCs according to a known 3D distribution
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

#ifndef __BCBENCHMARKMCMC3D__H
#define __BCBENCHMARKMCMC3D__H

#include <TF3.h>
#include <TH3F.h>

#include "BAT/BCModel.h"
#include "BAT/BCModelOutput.h"

//=============================================================================

class BCBenchmarkMCMC3D : public BCModel, public BCModelOutput
{
	public:

		BCBenchmarkMCMC3D(
				TF3* testFunction = NULL,
				const char* outputFile = "MCMCtest3d.root",
				const char* modelName = "BenchmarkMCMC3D");
		~BCBenchmarkMCMC3D();

		// inherited methods
		double LogAPrioriProbability(std::vector<double> parameters)
		{return 0;}

		double LogLikelihood(std::vector<double> parameters)
		{return log(fTestFunction->Eval(parameters[0],parameters[1],
				parameters[2]));}

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

		Double_t KolmogorovTest(TH1 *h1, TH1 *h2,
			Option_t *option);

		TH3F* GetTestFunctionHistogram();

		int fNbinx, fNbiny, fNbinz;
		double fXmin, fXmax, fYmin, fYmax, fZmin, fZmax;

		static const int fMaxChains = 5;
		static const int fMaxLags = 31; // lag = 1,2,...,30
		static const int fMax10thOfIters = 11; // iterations/10*(1,2,...10)

		TF3* fTestFunction;
		TF3* fFitFunction;

		TH3F* fHistXYZLags[fMaxChains][fMaxLags];
		TH3F* fHistXYZIter[fMaxChains][fMax10thOfIters];

		TH1F* fHChi2vsLags[fMaxChains];
		TH1F* fHChi2vsIter[fMaxChains];

		TH1F* fHKolmogorovProbVsLags[fMaxChains];
		TH1F* fHKolmogorovProbVsIter[fMaxChains];
};

#endif
