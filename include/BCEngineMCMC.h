/*!
 * \class BCEngineMCMC
 * \brief An engine class for Markov Chain Monte Carlo
 * \author Daniel Kollar
 * \author Kevin Kr&ouml;ninger
 * \version 1.0
 * \date 08.2008
 * \detail This class represents an engine class for Markov Chain
 * Monte Carlo (MCMC). One or more chains can be defined
 * simultaneously.
 */

/*
 * Copyright (C) 2008, Daniel Kollar and Kevin Kroeninger.
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

// ---------------------------------------------------------

#ifndef __BCENGINEMCMC__H
#define __BCENGINEMCMC__H

// ---------------------------------------------------------

#include <iostream>
#include <vector>

#include <math.h>

#include <TH1D.h>
#include <TH2D.h>
#include <TTree.h>
#include <TRandom3.h>
#include <TPrincipal.h>


// ---------------------------------------------------------

class BCEngineMCMC
{

	public:

		/** \name Constructors and destructors */
		/* @{ */

		/**
		 * Default constructor. */
		BCEngineMCMC();

		/**
		 * Constructor.
		 * @param n number of chains */
		BCEngineMCMC(int n);

		/**
		 * Default copy constructor. */
		BCEngineMCMC(const BCEngineMCMC & enginemcmc);

		/**
		 * Default destructor. */
		virtual ~BCEngineMCMC();

		/* @} */
		/** \name Assignment operators */
		/* @{ */

		/**
		 * Defaut assignment operator */
		BCEngineMCMC & operator = (const BCEngineMCMC & engineMCMC);

		/* @} */
		/** \name Getters */
		/* @{ */

		/*
		 * @return number of parameters of the Markov chain */
		int MCMCGetNParameters()
			{ return fMCMCNParameters; };

		/*
		 * @return number of Markov chains */
		int MCMCGetNChains()
			{ return fMCMCNChains; };

		/*
		 * @return number of iterations */
		std::vector <int> MCMCGetNIterations()
			{ return fMCMCNIterations; };

		/*
		 * @return pointer to the number of iterations */
		std::vector <int> * MCMCGetP2NIterations()
			{ return &fMCMCNIterations; };

		/*
		 * @return number of iterations needed for each chain to
		 * converge */
		std::vector <int> MCMCGetNIterationsConvergenceLocal()
			{ return fMCMCNIterationsConvergenceLocal; };

		/*
		 * @return number of iterations needed for all chains to
		 * converge simultaneously */
		int MCMCGetNIterationsConvergenceGlobal()
			{ return fMCMCNIterationsConvergenceGlobal; };

		/*
		 * @return flag if converged or not */
		bool MCMCGetFlagConvergenceGlobal()
			{ return fMCMCFlagConvergenceGlobal; };

		/*
		 * @return maximum number of iterations for a Markov chain */
		int MCMCGetNIterationsMax()
			{ return fMCMCNIterationsMax; };

		/*
		 * @return number of iterations needed for burn-in. These
		 * iterations are not included in fMCMCNIterations */
		int MCMCGetNIterationsBurnIn()
			{ return fMCMCNIterationsBurnIn; };

		/*
		 * @return number of iterations needed for PCA. These
		 * iterations are not included in fMCMCNIterations */
		int MCMCGetNIterationsPCA()
			{ return fMCMCNIterationsPCA; };

		/*
		 * @returns number of accepted trials for each chain */
		std::vector <int> MCMCGetNTrialsTrue()
			{ return fMCMCNTrialsTrue; };

		/*
		 * @returns number of not-accepted trials for each chain */
		std::vector <int> MCMCGetNTrialsFalse()
			{ return fMCMCNTrialsFalse; };

		/*
		 * @return mean value of the probability for each chain up to
		 * the current iteration  */
		std::vector <double> MCMCGetMean()
			{ return fMCMCMean; };

		/*
		 * @return mean value of the probability for each chain up to
		 * the current iteration */
		std::vector <double> MCMCGetVariance()
			{ return fMCMCVariance; };

		/*
		 * @return flag to automatically calculate the number of iterations
		 * of a Markov chain */
		bool MCMCGetFlagIterationsAuto()
			{ return fMCMCFlagIterationsAuto; };

		/*
		 * @return scale factor for the width of the trial function */
		double MCMCGetTrialFunctionScale()
			{ return fMCMCTrialFunctionScale; };

		/*
		 * @return scale factor for all parameters and chains */
		std::vector <double> MCMCGetTrialFunctionScaleFactor()
			{ return fMCMCTrialFunctionScaleFactor; };

		/*
		 * @return scale factor for all parameters and achain.
		 * @param ichain chain index */
		std::vector <double> MCMCGetTrialFunctionScaleFactor(int ichain);

		/*
		 * @return scale factor for a parameter and a chain.
		 * @param ichain chain index
		 * @param ipar parameter index */
		double MCMCGetTrialFunctionScaleFactor(int ichain, int ipar);

		/*
		 * @return current point of each Markov chain */
		std::vector <double> MCMCGetx()
			{ return fMCMCx; };

		/*
		 * @return pointer of each Markov chain */
		std::vector <double> * MCMCGetP2x()
			{ return &fMCMCx; };

		/*
		 * @param ichain index of the Markov chain
		 * @return current point of the Markov chain */
		std::vector <double> MCMCGetx(int ichain);

		/*
		 * @param ichain chain index
		 * @param ipar parameter index
		 * @return parameter of the Markov chain */
		double MCMCGetx(int ichain, int ipar);

		/*
		 * @return log of the probability of the current points of each Markov chain */
		std::vector <double> MCMCGetLogProbx()
			{ return fMCMCLogProbx; };

		/*
		 * @return log of the probability of the current points of the Markov chain.
		 * @param ichain chain index */
		double MCMCGetLogProbx(int ichain);

		/*
		 * @return pointer to the log of the probability of the current points of each Markov chain */
		std::vector <double> * MCMCGetP2LogProbx()
			{ return &fMCMCLogProbx; };

		/*
		 * @return maximum points of each Markov chain */
		std::vector <double> MCMCGetMaximumPoints()
			{ return fMCMCMaximumPoints; };

		/*
		 * @return maximum point of  Markov chain
		 * @param i The index of the Markov chain */
		std::vector <double> MCMCGetMaximumPoint(int i);

		/*
		 * @return maximum (log) probability of each Markov chain */
		std::vector <double> MCMCGetMaximumLogProb()
			{ return fMCMCMaximumLogProb; };

		/*
		 * @return control plots */
		//	TH1D * MCMCGetH1RValue()
		//	  { return fMCMCH1RValue; };

		//	TH1D * MCMCGetH1Efficiency()
		//	  { return fMCMCH1Efficiency; };

		/*
		 * @return flag which defined initial position */
		int MCMCGetFlagInitialPosition()
			{ return fMCMCFlagInitialPosition; };

		/*
		 * @return R-value criterion */
		double MCMCGetRValueCriterion()
			{ return fMCMCRValueCriterion; };

		/*
		 * @return R-value criterion for parameters */
		double MCMCGetRValueParametersCriterion()
			{ return fMCMCRValueParametersCriterion; };

		/*
		 * @return R-value */
		double MCMCGetRValue()
			{ return fMCMCRValue; };

		/*
		 * @return R-value for a parameter
		 * @param i parameter index */
		double MCMCGetRValueParameters(int i)
			{ return fMCMCRValueParameters.at(i); };

		/*
		 * @return the relative precision for the estimate of the mode */
//		double MCMCGetPrecisionMode()
//			{ return fMCMCRelativePrecisionMode; };

		/*
		 * @return the flag for the use of PCA */
		bool MCMCGetFlagPCA()
			{ return fMCMCFlagPCA; };

		/*
		 * Rtrieve the tree containing the Markov chain.
		 * @param i index of the Markov chain
		 * @return pointer to the tree */
		TTree * MCMCGetMarkovChainTree(int i)
			{ return fMCMCTrees.at(i); };

		/*
		 * Retrieve a histogram of the 1D marginalized distribution of a single parameter.
		 * @param i index of the parameter
		 * @return pointer to the histogram */
		TH1D * MCMCGetH1Marginalized(int i)
			{ return fMCMCH1Marginalized[i]; };

		/*
		 * Retrieve a histogram of the 2D marginalized distribution for two parameters.
		 * @param i index of the first parameter
		 * @param j index of the second parameter
		 * @return pointer to the histogram */
		TH2D * MCMCGetH2Marginalized(int i, int j);

		/* @} */
		/** \name Setters */
		/* @{ */

		/*
		 * Sets the scale factor for the width of the trial function. */
		void MCMCSetTrialFunctionScale(double scale)
			{ fMCMCTrialFunctionScale = scale; };

		void MCMCSetTrialFunctionScaleFactor(std::vector <double> scale)
			{ fMCMCTrialFunctionScaleFactorStart = scale; };

		/*
		 * Sets the number of parameters of the Markov chain. */
		void MCMCSetNParameters(int n);

		/*
		 * Sets the number of Markov chains which are run in parallel. */
		void MCMCSetNChains(int n);

		/*
		 * Sets the maximum number of iterations. */
		void MCMCSetNIterationsMax(int n)
			{ fMCMCNIterationsMax = n; };

		/*
		 * Sets the number of iterations. */
		void MCMCSetNIterationsRun(int n)
			{ fMCMCNIterationsRun = n; };

		/*
		 * Sets the number of iterations needed for burn-in. */
		void MCMCSetNIterationsBurnIn(int n)
			{ fMCMCNIterationsBurnIn = n; };

		/*
		 * Sets the number of iterations needed for PCA. */
		void MCMCSetNIterationsPCA(int n)
			{ fMCMCNIterationsPCA = n; };

		void MCMCSetNIterationsUpdate(int n)
		{ fMCMCNIterationsUpdate = n; };

		/*
		 * Sets flag to automatically calculate the number of iterations of
		 * a Markov chain. */
		void MCMCSetIterationsAuto(bool flag)
			{ fMCMCFlagIterationsAuto = flag; };

		/*
		 * Sets the minimum efficiency required for a chain. */
		void MCMCSetMinimumEfficiency(double efficiency)
			{ fMCMCEfficiencyMin = efficiency; };

		/*
		 * Sets the maximum efficiency required for a chain. */
		void MCMCSetMaximumEfficiency(double efficiency)
			{ fMCMCEfficiencyMax = efficiency; };

		/*
		 * Sets flag to write Markov chains to file. */
		void MCMCSetWriteChainToFile(bool flag)
			{ fMCMCFlagWriteChainToFile = flag; };

		/*
		 * Sets the initial position for a chain.
		 * @param chain chain index.
		 * @param x0 intial position
		 * @see MCMCSetIntitialPositions. */
		void MCMCSetInitialPosition(std::vector<double> x0);
		void MCMCSetInitialPosition(int chain, std::vector<double> x0);

		/*
		 * Sets the initial positions for all chains.
		 * @param x0s initial positions for all chains */
		void MCMCSetInitialPositions(std::vector<double> x0s);

		/*
		 * Sets flag which defined initial position.
		 */
		void MCMCSetFlagInitialPosition(int flag)
			{ fMCMCFlagInitialPosition = flag; };

		/*
		 * Sets the flag which controls the sequence parameters during the running
		 * of the MCMC.  */
		void MCMCSetFlagOrderParameters(bool flag)
		{ fMCMCFlagOrderParameters = flag; };

		/*
		 * Sets the R-value criterion for convergence of all chains. */
		void MCMCSetRValueCriterion(double r)
			{ fMCMCRValueCriterion = r; };

		/*
		 * Sets the parameter R-value criterion for convergence of all chains */
		void MCMCSetRValueParametersCriterion(double r)
			{ fMCMCRValueParametersCriterion = r; };

		/*
		 * Sets the relative precision for the estimate of the mode. */
//		void MCMCSetPrecisionMode(double precision)
//			{ fMCMCRelativePrecisionMode = precision; };

		/*
		 * Sets the flag to either perform a pre-run with PCA or not. */
		void MCMCSetFlagPCA(bool flag)
			{ fMCMCFlagPCA = flag; };

		/*
		 * Sets the tree containing the Markov chains. */
		void MCMCSetMarkovChainTrees(std::vector <TTree *> trees);

		/*
		 * Set a flag to control if during the PCA the least eigenvectors
		 * should be ignored or not. */
		void MCMCSetFlagPCATruncate(bool flag)
			{ fMCMCFlagPCATruncate = flag; };

		/*
		 * Sets the minimum ratio of an eigenvalue to the largest eigenvalue
		 * below which it is ignored if fMCMCFlagPCATruncate is true. */
		void MCMCSetPCAMinimumRatio(double ratio)
			{ fMCMCPCAMinimumRatio = ratio; };

		/* @} */
		/** \name Miscellaneous methods */
		/* @{ */

		/*
		 * Adds a parameter.
		 * @param min minimum value of the parameter
		 * @param max maximum value of the parameter
		 * @return number of parameters after adding */
		int MCMCAddParameter(double min, double max);

		/*
		 * Random walk trial function. The function is symmetric and
		 * used for the Metropolis algorithm.
		 * @param x point with the dimension fMCMCNParameters */
		void MCMCTrialFunction(std::vector <double> &x);
		void MCMCTrialFunctionSingle(int ichain, int iparameter, std::vector <double> &x);

		/*
		 * Independent chain trial function. The function does not
		 * have to be symmetric and is used for the
		 * Metropolis-Hastings algorithm.
		 * @param x point with the dimension fMCMCNParameters
		 * @return transition probability */
		double MCMCTrialFunctionIndependent(std::vector <double> &xnew, std::vector <double> &xold, bool newpoint);

		/*
		 * Trial function.
		 * @param x point with the dimension fMCMCNParameters */
		void MCMCTrialFunctionAuto(std::vector <double> &x);

		/*
		 * Trial function for the MCMC relative to the old point. No PCA.
		 * @param pxold pointer to the old point. The length of the vector equals to fMCMCNParameters.
		 * @param pxnew pointer to the new point. The length of the vector equals to fMCMCNParameters.
		 * @param flag_compute flag which indicates whether to compute a new point (true) or to just pass the value of the function (flase)
		 * @return value of the trial function */
		double MCMCTrialFunctionRelativeNoPCA(std::vector <double> * xold, std::vector<double> * xnew, bool flag_compute);

		/*
		 * not documented !!! */
		void MCMCGetProposalPoint(int chain, std::vector <double> xnew, std::vector <double> xold);

		/*
		 * Returns a trial point for the Metropolis algorithm.
		 * @param chain chain index
		 * @param x proposal point
		 * @param pca bool whether to use PCA or not
		 * @return flag indicating whether the new point lies within the allowed range */
		bool MCMCGetProposalPointMetropolis(int chain, std::vector <double> &x, bool pca);

		/*
		 * Returns a trial point for the Metropolis algorithm.
		 * @param chain chain index
		 * @param x proposal point
		 * @param pca bool whether to use PCA or not
		 * @return flag indicating whether the new point lies within the allowed range */
		bool MCMCGetProposalPointMetropolis(int chain, int parameter, std::vector <double> &x, bool pca);

		/*
		 * Returns a trial point for the Metropolis algorithm.
		 * @param chain chain index
		 * @param x proposal point
		 * @param pca bool whether to use PCA or not
		 * @return flag indicating whether the new point lies within the allowed range */
		bool MCMCGetProposalPointMetropolisHastings(int chain, std::vector <double> &xnew, std::vector <double> &xold);

		/*
		 * This method samples uniformly over the allowed parameter space.
		 * @return new point for the PCA run */
		void MCMCGetNewPointPCA();

		/*
		 * Generates a new point using the Metropolis algorithm.
		 * @param chain chain index
		 * @param pca bool whether to use PCA or not */
		bool MCMCGetNewPointMetropolis(int chain = 0, bool pca = false);
		bool MCMCGetNewPointMetropolis(int chain = 0, int parameter = 0, bool pca = false);

		/*
		 * Generates a new point using the Metropolis algorithm.
		 * @param chain chain index */
		bool MCMCGetNewPointMetropolisHastings(int chain = 0);

		/*
		 * Generates a new point using simulated annealing.
		 * @param chain chain index
		 * @param pca bool whether to use PCA or not */
		bool MCMCGetNewPointSimulatedAnnealing(int chain = 0, bool pca = false);

		/*
		 * Calculates the temperature accoring to the annealing schedule.
		 * @param chain chain index
		 * @return temperature */
		double MCMCAnnealingSchedule(int chain);

		/*
		 * Updates statistics, plots, etc. */
		void MCMCUpdateStatistics();
		void MCMCUpdateStatisticsModeConvergence();
		void MCMCUpdateStatisticsCheckMaximum();
		void MCMCUpdateStatisticsFillHistograms();
		void MCMCUpdateStatisticsTestConvergenceAllChains();
		void MCMCUpdateStatisticsWriteChains();

		/*
		 * Needs to be overloaded in the derived class.
		 * @return natural logarithm of the function to map with MCMC */
		virtual double LogEval(std::vector <double> parameters);

		/*
		 * Perform a run for the PCA */
//		void MCMCPCARun();

		/*
		 * Runs Metropolis algorithm. */
		int MCMCMetropolis();

		/*
		 * Runs a pre run for the Metropolis algorithm. */
		int MCMCMetropolisPreRun();

		/*
		 * Runs Metropolis-Hastings algorithm. */
		int MCMCMetropolisHastings();

		/*
		 * Runs simulated annealing algorithm. */
		int MCMCSimulatedAnnealing();

		/*
		 * Performs an initial run. */
		void MCMCInitialRun();

		/*
		 * Resets the run statistics. */
		void MCMCResetRunStatistics();

		/*
		 * Initializes Markov chains. */
		void MCMCInitializeMarkovChains();

		/*
		 * Initializes the engine. */
		int MCMCInitialize();

		/*
		 * Interface allowing to execute arbitrary code for each iteration of the MCMC.
		 * This method needs to be overloaded in the derived class. */
		virtual void MCMCIterationInterface()
			{};

		/*
		 * Sets the histogram with 1D marginalized distributions for parameter.
		 * @param i index of the parameter
		 * @param h pointer to an existing histogram */
		int SetMarginalized(int index, TH1D * h);

		/*
		 * Sets the histogram with 2D marginalized distributions for two parameters.
		 * @param index1 index of the first parameter
		 * @param index2 index of the second parameter
		 * @param h pointer to an existing histogram */
		int SetMarginalized(int index1, int index2, TH2D * h);

		/**
		 * Set the default values for the MCMC chain. */
		void MCMCSetValuesDefault();

		/**
		 * Set the values for a quick MCMC run. */
		void MCMCSetValuesQuick();

		/**
		 * Set the values for a detailed MCMC run. */
		void MCMCSetValuesDetail();

		/* @} */

	private:

		/*
		 * Copies this BCEngineMCMC into another one. */
		void Copy(BCEngineMCMC & enginemcmc) const;

		/*
		 * Defines a type of a pointer to a member function. */
		typedef bool (BCEngineMCMC::*MCMCPointerToGetProposalPoint) (int chain, std::vector <double> xnew, std::vector <double> xold) const;

		/*
		 * Pointer to a member function */
		MCMCPointerToGetProposalPoint fMCMCPointerToGetProposalPoint;

	protected:

		/*
		 * Number of parameters */
		int fMCMCNParameters;

		/*
		 * Parameter boundaries */
		std::vector <double> fMCMCBoundaryMin;
		std::vector <double> fMCMCBoundaryMax;

		/*
		 * Number of Markov chains ran in parallel */
		int fMCMCNChains;

		/*
		 * Number of total iterations of the Markov chains. The length of
		 * the vector is equal to fMCMCNChains. */
		std::vector<int> fMCMCNIterations;

		/*
		 * Number of iterations for updating scale factors */
		int fMCMCNIterationsUpdate;

		/*
		 * Number of iterations needed for each chain to convergence. The
		 * size of the vector is equal to fMCMCNChains. */
		std::vector<int> fMCMCNIterationsConvergenceLocal;

		/*
		 * Number of iterations needed for all chains to convergence simulaneously */
		int fMCMCNIterationsConvergenceGlobal;

		/*
		 * Flag for convergence */
		bool fMCMCFlagConvergenceGlobal;

		/*
		 * Maximum number of iterations for a Markov chain prerun */
		int fMCMCNIterationsMax;

		/*
		 * Number of iterations for a Markov chain run */
		int fMCMCNIterationsRun;

		/*
		 * Minimum number of iterations for the pre-run */
		int fMCMCNIterationsPreRunMin;

		/*
		 * Number of iterations for burn-in. These iterations are not
		 * included in fMCMCNIterations. */
		int fMCMCNIterationsBurnIn;

		/*
		 * Number of iterations for PCA. These iterations are not included
		 * in FMCMCNIterations. */
		int fMCMCNIterationsPCA;

		/*
		 * Mean values of the data in the PCA coordinate system */
		std::vector <double> fMCMCPCAMean;

		/*
		 * Variance of the data in the PCA coordinate system */
		std::vector <double> fMCMCPCAVariance;

		/*
		 * Flag to control if during the PCA the least eigenvectors should
		 * be ignored or not */
		bool fMCMCFlagPCATruncate;

		/*
		 * Minimum ratio of an eigenvalue to the largest eigenvalue below
		 * which it is ignored if fMCMCFlagPCATruncate is true */
		double fMCMCPCAMinimumRatio;

		/*
		 * If the least eigenvectors are ignored this is the number of
		 * dimensions remaining */
		int fMCMCNDimensionsPCA;

		/*
		 * Number of accepted trials and not accepted trials for each
		 * chain. The length of the vectors is equal to fMCMCNChains *
		 * fMCMCNParameters. For each chain these numbers add up to
		 * fMCMCNIterations. */
		std::vector<int> fMCMCNTrialsTrue;
		std::vector<int> fMCMCNTrialsFalse;

		/*
		 * The mean and variance of all values of each Markov chain. The
		 * length of the vectors is equal to fMCMCNChains. */
		std::vector <double> fMCMCMean;
		std::vector <double> fMCMCVariance;

		/*
		 * The sum and sum squared of all values of each Markov chain. These
		 * are used to calculate the mean and variance of each chain. The
		 * length of the vectors is equal to fMCMCNChains. */
		std::vector <double> fMCMCSum;
		std::vector <double> fMCMCSum2;

		/*
		 * The mean and variance of all parameters of each Markov chain. The
		 * length of the vectors is equal to fMCMCNChains * fMCMCNParameters. */
		std::vector <double> fMCMCMeanx;
		std::vector <double> fMCMCVariancex;

		/*
		 * Flag to automatically calculate the number of iterations of a
		 * Markov chain */
		bool fMCMCFlagIterationsAuto;

		/*
		 * Flag to write Markov chains to file */
		bool fMCMCFlagWriteChainToFile;

		/*
		 * Scales the width of the trial functions by a global factor */
		double fMCMCTrialFunctionScale;

		/*
		 * Scales the width of the trial functions by a scale factor
		 * for each parameter and chain */
		std::vector <double> fMCMCTrialFunctionScaleFactor;
		std::vector <double> fMCMCTrialFunctionScaleFactorStart;

		/*
		 * Defines if a prerun has been performed or not */
		bool fMCMCFlagPreRun;

		/*
		 * The intial position of each Markov chain. The length of the
		 * vectors is equal to fMCMCNChains * fMCMCNParameters. First, the
		 * values of the first Markov chain are saved, then those of the
		 * second and so on */
		std::vector <double> fMCMCInitialPosition;

		/*
		 * Minimum efficiency for MCMC */
		double fMCMCEfficiencyMin;

		/*
		 * Maximum efficiency for MCMC */
		double fMCMCEfficiencyMax;

		/*
		 * Variable which defines the initial position. 0 (default) center
		 * of the allowed region, (1) random initial position (2)
		 * pre-defined intial position. */
		int fMCMCFlagInitialPosition;

		/*
		 * Flag which controls the sequence parameters during the running 
		 * of the MCMC. 
		 */ 
		bool fMCMCFlagOrderParameters; 

		/*
		 * The current points of each Markov chain. The length of the
		 * vectors is equal to fMCMCNChains * fMCMCNParameters. First, the
		 * values of the first Markov chain are saved, then those of the
		 * second and so on. */
		std::vector <double> fMCMCx;

		/*
		 * A temporary vector for a single Markov chain */
		std::vector <double> fMCMCxLocal;

		/*
		 * The log of the probability of the current points of each Markov
		 * chain. The length of the vectors is fMCMCNChains. */
		std::vector<double> fMCMCLogProbx;

		/*
		 * The maximum points of each Markov chain. The length of the vector
		 * is fMCMCNChains * fMCMCNParameters. First, the values of the
		 * first Markov chain are saved, then those of the second and so on. */
		std::vector <double> fMCMCMaximumPoints;

		/*
		 * The maximum (log) probability of each Markov chain. The length of
		 * the vector is fMCMCNChains. */
		std::vector <double> fMCMCMaximumLogProb;

		/*
		 * The R-value criterion for convergence of log-likelihood*/
		double fMCMCRValueCriterion;

		/*
		 * The R-value criterion for convergence of parameters */
		double fMCMCRValueParametersCriterion;

		/*
		 * The R-value at which the chains did converge */
		double fMCMCRValue;

		std::vector <double> fMCMCRValueParameters;

		/*
		 * The relative precision for the convergence of the modes of
		 * several chains */
//		double fMCMCRelativePrecisionMode;

//		std::vector <double> fMCMCRelativePrecisionModeValues;
		std::vector <double> fMCMCNumericalPrecisionModeValues;

		/*
		 * The starting temperature for simulated annealing */
		double fMCMCSimulatedAnnealingT0;

		/*
		 * Random number generator */
		TRandom3 * fMCMCRandom;

		/*
		 * PCA */
		TPrincipal * fMCMCPCA;

		/*
		 * Flag to perform PCA or not */
		bool fMCMCFlagPCA;

		/*
		 * Number of bins per dimension for the marginalized distributions. */
		std::vector<int> fMCMCH1NBins;

		/*
		 * An array of marginalized distributions */
		std::vector <TH1D *> fMCMCH1Marginalized;
		std::vector <TH2D *> fMCMCH2Marginalized;

		/*
		 * Control plots */
//		TH1D * fMCMCH1RValue;     // R-value
//		TH1D * fMCMCH1Efficiency; // Efficiency of the Markov chain

		/*
		 * The trees containing the Markov chains. The length of the vector
		 * is fMCMCNChains. */
		std::vector<TTree *> fMCMCTrees;
};

// ---------------------------------------------------------

#endif
