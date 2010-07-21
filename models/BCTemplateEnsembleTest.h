#ifndef __BCTEMPLATEENSEMBLETEST__H
#define __BCTEMPLATEENSEMBLETEST__H

/*!
 * \class BCTemplateEnsembleTest
 * This class can be used for ensemble tests using the StackTool. The
 * fitting can be done with Minuit or with Markov Chains.
 *
 * \brief A class for doing ensemble tests.
 * \author Andrea Knue
 * \author Daniel Kollar
 * \author Kevin Kroeninger
 * \date 10.04.2010
 */

/*
 * Copyright (C) 2008, 2009, 2010, Daniel Kollar and Kevin Kroeninger.
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

// --------------------------------------------------------------------------------------------

#include <TH1D.h>

class TFile;
class TTree;
class TRandom3;

class BCTemplateFitter;

// --------------------------------------------------------------------------------------------

class BCTemplateEnsembleTest
{
	public:

	  /**
	   * The constructor.
	   */
	  BCTemplateEnsembleTest();

	  /**
	   * The destructor.
	   */
	  ~BCTemplateEnsembleTest();

	  /**
	   * Set the template used to generate the ensembles.
	   */
	  int SetEnsembleTemplate(TH1D hist);

		/**
		 * Set the BCTemplateFitter used to analyze the ensembles.
		 */
		void SetTemplateFitter(BCTemplateFitter * model)
			{ fTemplateFitter = model; }

	  /**
	   * A function to perform an ensemble test for each data set in the container.
	   */
	  int PerformEnsembleTest();

	  /**
	   * A function to define the number of events per ensemble.
	   */
	  void SetEnsembleExpectation(double expectation)
			{ fEnsembleExpectation = expectation; }

	  /**
	   * A function to define the number of ensembles per data set.
	   */
	  void SetNEnsembles(int n)
			{ fNEnsembles = n; }

	  /**
	   * A function to set the MCMC flag.
	   */
	  void SetFlagMCMC(bool flag)
	  { fFlagMCMC = false; } // debugKK

		/**
		 * Write tree to file
		 */
		int Write(const char * filename);

		/**
		 * Prepare tree.
		 */
		int PrepareTree();

	private:

	  /**
	   * Create a new ensemble.
		 * @return A histogram with the new ensemble.
	   */
	  TH1D * BuildEnsemble();

	protected:

	  /**
	   * The template used for the generation of ensembles.
	   */
		TH1D fEnsembleTemplate;

		/**
		 * The stack model used to analyze the ensembles.
		 */
		BCTemplateFitter * fTemplateFitter;

		/**
		 * Output file.
		 */
		TFile * fFile;

		/**
		 * Output tree.
		 */
		TTree * fTree;

		/**
		 * A counter for the number of ensembles.
		 */
		int fEnsembleCounter;

	  /**
	   * Exepectation value
	   */
	  double fEnsembleExpectation;

	  /**
	   * Number of ensembles per data set.
	   */
	  int fNEnsembles;

	  /**
	   * A flag to turn the Markov Chains on.
	   */
	  bool fFlagMCMC;

		/**
		 * The random number generator.
		 */
		TRandom3 * fRandom;

		/**
		 * Tree variable: global mode
		 */
		std::vector<double> fOutParModeGlobal;

		/**
		 * Tree variable: positive uncertainty on global mode
		 */
		std::vector<double> fOutParErrorUpGlobal;

		/**
		 * Tree variable: negative uncertainty on global mode
		 */
		std::vector<double> fOutParErrorDownGlobal;

		/**
		 * Tree variable: marginalized mode
		 */
		std::vector<double> fOutParModeMarg;

		/**
		 * Tree variable: median
		 */
		std::vector<double> fOutParMedianMarg;

		/**
		 * Tree variable: mean
		 */
		std::vector<double> fOutParMeanMarg;

		/**
		 * Tree variable: rms
		 */
		std::vector<double> fOutParRMSMarg;

		/**
		 * Tree variable: 16% quantile
		 */
		std::vector<double> fOutParErrorUpMarg;

		/**
		 * Tree variable: 84% quantile
		 */
		std::vector<double> fOutParErrorDownMarg;

		/**
		 * Tree variable: 5% quantile
		 */
		std::vector<double> fOutParQuantile5Marg;

		/**
		 * Tree variable: 10% quantile
		 */
		std::vector<double> fOutParQuantile10Marg;

		/**
		 * Tree variable: 90% quantile
		 */
		std::vector<double> fOutParQuantile90Marg;

		/**
		 * Tree variable: 95% quantile
		 */
		std::vector<double> fOutParQuantile95Marg;

		/**
		 * Tree variable: marginalized mode (ratio)
		 */
		std::vector<double> fOutRatioModeMarg;

		/**
		 * Tree variable: median (ratio)
		 */
		std::vector<double> fOutRatioMedianMarg;

		/**
		 * Tree variable: mean (ratio)
		 */
		std::vector<double> fOutRatioMeanMarg;

		/**
		 * Tree variable: rms (ratio)
		 */
		std::vector<double> fOutRatioRMSMarg;

		/**
		 * Tree variable: 16% quantile (ratio)
		 */
		std::vector<double> fOutRatioErrorUpMarg;

		/**
		 * Tree variable: 84% quantile (ratio)
		 */
		std::vector<double> fOutRatioErrorDownMarg;

		/**
		 * Tree variable: 5% quantile (ratio)
		 */
		std::vector<double> fOutRatioQuantile5Marg;

		/**
		 * Tree variable: 10% quantile (ratio)
		 */
		std::vector<double> fOutRatioQuantile10Marg;

		/**
		 * Tree variable: 90% quantile (ratio)
		 */
		std::vector<double> fOutRatioQuantile90Marg;

		/**
		 * Tree variable: 95% quantile (ratio)
		 */
		std::vector<double> fOutRatioQuantile95Marg;

		/**
		 * Tree variable: chi2
		 */
		double fOutChi2;

		/**
		 * Tree variable: ndf
		 */
		int fOutNDF;

		/**
		 * Tree variable: chi2-probability
		 */
		double fOutChi2Prob;

		/**
		 * Tree variable: KL probability
		 */
		double fOutKSProb;

		/**
		 * Tree variable: p-value
		 */
		double fOutPValue;

		/**
		 * Tree variable: number of events in the data
		 */
		int fOutNEvents;
};

#endif
