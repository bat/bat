/**  
 * \class BCEngineMCMC
 * \brief An eninge class for Markov Chain Monte Carlo 
 * \author D. Kollar  
 * \author K. Kr&ouml;ninger  
 * \version 1.0  
 * \date 16.12.2007  
 *  
 * This class defines an engine for Markov Chain Monte Carlo
 * (MCMC). One or more chains can be handled. 
 *  
 * Copyright (C) 2007, D. Kollar, K. Kr&ouml;ninger  
 */  

// --------------------------------------------------------- 

#ifndef __BCENGINEMCMC__H
#define __BCENGINEMCMC__H

// ---------------------------------------------------------

#include <iostream.h>
#include <math.h>

#include <TH1D.h> 
#include <TH2D.h> 
#include <TTree.h> 
#include <TRandom3.h> 
#include <TPrincipal.h> 

#include <vector.h>

// --------------------------------------------------------- 

class BCEngineMCMC
{

 public:   

	/** \name Constructors and destructors */ 
	/* @{ */

	/** 
	 * The default constructor. 
	 * @param n Number of chains. 
	 */ 
	BCEngineMCMC(); 

	/** 
	 * A constructor 
	 */ 
	BCEngineMCMC(int n); 

	/** 
	 * The default copy constructor. 
	 */ 
	BCEngineMCMC(const BCEngineMCMC & enginemcmc); 

	/** 
	 * The default destructor. 
	 */ 
	virtual ~BCEngineMCMC(); 

	/* @} */ 

	/** \name Assignment operators */ 
	/* @{ */ 

	/**
	 * The defaut assignment operator 
	 */ 
	BCEngineMCMC & operator = (const BCEngineMCMC & engineMCMC); 

	/* @} */ 

	/** \name Getters */ 
	/* @{ */ 

	/*
	 * Returns the number of parameters of the Markov chain. 
	 */ 
	int MCMCGetNParameters()
	{ return fMCMCNParameters; }; 

	/*
	 * Returns the number of Markov chains. 
	 */ 
	int MCMCGetNChains()
	{ return fMCMCNChains; }; 

	/*
	 * Returns the number of iterations. 
	 */
	std::vector <int> MCMCGetNIterations()
		{ return fMCMCNIterations; }; 

	/*
	 * Returns a pointer to the number of iterations. 
	 */
	std::vector <int> * MCMCGetP2NIterations()
		{ return &fMCMCNIterations; }; 

	/*
	 * Returns the number of iterations needed for each chain to
	 * converge. 
	 */ 
	std::vector <int> MCMCGetNIterationsConvergenceLocal()
		{ return fMCMCNIterationsConvergenceLocal; }; 

	/*
	 * Returns the number of iterations needed for all chains to
	 * converge simultaneously.
	 */ 
	int MCMCGetNIterationsConvergenceGlobal()
	{ return fMCMCNIterationsConvergenceGlobal; }; 

	/*
	 * Returns a flag if converged or not. 
	 */ 
	bool MCMCGetFlagConvergenceGlobal()
	{ return fMCMCFlagConvergenceGlobal; }; 

	/*
	 * Returns the maximum number of iterations for a Markov chain. 
	 */ 
	int MCMCGetNIterationsMax() 
	{ return fMCMCNIterationsMax; }; 

	/*
	 * Returns the number of iterations needed for burn-in. These
	 * iterations are not included in fMCMCNIterations. 
	 */ 
	int MCMCGetNIterationsBurnIn()
	{ return fMCMCNIterationsBurnIn; }; 

	/*
	 * Returns the number of iterations needed for PCA. These
	 * iterations are not included in fMCMCNIterations. 
	 */ 
	int MCMCGetNIterationsPCA()
	{ return fMCMCNIterationsPCA; }; 

	/*
	 * Returns the number of accepted trials for each chain.
	 */ 
	std::vector <int> MCMCGetNTrialsTrue()
		{ return fMCMCNTrialsTrue; }; 

	/*
	 * Returns the number of not-accepted trials for each chain.
	 */ 
	std::vector <int> MCMCGetNTrialsFalse()
		{ return fMCMCNTrialsFalse; }; 

	/* 
	 * Returns the mean value of the probability for each chain up to
	 * the current iteration. 
	 */ 
	std::vector <double> MCMCGetMean() 
		{ return fMCMCMean; }; 

	/* 
	 * Returns the mean value of the probability for each chain up to
	 * the current iteration. 
	 */ 
	std::vector <double> MCMCGetVariance()
		{ return fMCMCVariance; }; 

	/*
	 * Returns flag to automatically calculate the number of iterations
	 * of a Markov chain.
	 */ 
	bool MCMCGetFlagIterationsAuto()
	{ return fMCMCFlagIterationsAuto; }; 

	/*
	 * Returns the scale factor for the width of the trial function.
	 */ 
	double MCMCGetTrialFunctionScale()
	{ return fMCMCTrialFunctionScale; }; 

	/*
	 * Returns the scale factor for all parameters and chains. 
	 */ 
	std::vector <double> MCMCGetTrialFunctionScaleFactor()
	{ return fMCMCTrialFunctionScaleFactor; }; 

	/*
	 * Returns the scale factor for all parameters and the ith
	 * chain.
	 * @param i The chain. 
	 */ 
	std::vector <double> MCMCGetTrialFunctionScaleFactor(int i); 

	/*
	 * Returns the scale factor for the jth parameter and the ith
	 * chain.  
	 * @param i The chain.
	 * @param j The parameter index.
	 */ 
	double MCMCGetTrialFunctionScaleFactor(int i, int j); 

	/*
	 * Returns the current point of each Markov chain. 
	 */ 
	std::vector <double> MCMCGetx()
		{ return fMCMCx; }; 

	/* 
	 * Returns the pointer of each Markov chain.
	 */ 
	std::vector <double> * MCMCGetP2x() 
		{ return &fMCMCx; }; 

	/* 
	 * Returns the current point of the ith Markov chain.
	 * @param i The index of the Markov chain. 
	 */ 
	std::vector <double> MCMCGetx(int i); 

	/* 
	 * Returns the jth parameter of the ith Markov chain.
	 * @param i The index of the Markov chain. 
	 * @param j The parameter index
	 */ 
	double MCMCGetx(int i, int j); 

	/* 
	 * Returns the log of the probability of the current points of each
	 * Markov chain. 
	 */ 
	std::vector <double> MCMCGetLogProbx()
		{ return fMCMCLogProbx; }; 

	/* 
	 * Returns the log of the probability of the current points of the
	 * ith Markov chain.
	 * @param i The index of the Markov chain
	 */ 
	double MCMCGetLogProbx(int i);

	/* 
	 * Returns the pointer to the log of the probability of the current
	 * points of each Markov chain.
	 */ 
	std::vector <double> * MCMCGetP2LogProbx()
		{ return &fMCMCLogProbx; }; 

	/*
	 * Returns the maximum points of each Markov chain. 
	 */ 
	std::vector <double> MCMCGetMaximumPoints()
		{ return fMCMCMaximumPoints; }; 

	/*
	 * Returns the maximum point of each ith Markov chain. 
	 * @param i The index of the Markov chain. 
	 */ 
	std::vector <double> MCMCGetMaximumPoint(int i); 

	/*
	 * Returns the maximum (log) probability of each Markov chain.
	 */ 
	std::vector <double> MCMCGetMaximumLogProb()
		{ return fMCMCMaximumLogProb; }; 

	/*
	 * Returns the control plots 
	 */ 
	TH1D * MCMCGetH1RValue()
	  { return fMCMCH1RValue; }; 

	TH1D * MCMCGetH1Efficiency()
	  { return fMCMCH1Efficiency; }; 

	/*
	 * Returns flag which defined initial position. 
	 */ 
	int MCMCGetFlagInitialPosition() 
	{ return fMCMCFlagInitialPosition; }; 

	/*
	 * Returns the r-value criterion 
	 */ 
	double MCMCGetRValueCriterion()
	{ return fMCMCRValueCriterion; }; 

	/*
	 * Returns the r-value criterion 
	 */ 
	double MCMCGetRValue()
	{ return fMCMCRValue; }; 

	/*
	 * Returns the flag for the use of PCA 
	 */ 
	bool MCMCGetFlagPCA()
	{ return fMCMCFlagPCA; }; 

	/*
	 * Returns the tree containing the ith Markov chain. 
	 * @param i The index of the Markov chain.
	 * @return A pointer to the tree. 
	 */ 
	TTree * MCMCGetMarkovChainTree(int i)
		{ return fMCMCTrees.at(i); }; 

	/*
	 * Returns a histogram of the 1-d marginalized distribution of the ith
	 * parameter. 
	 * @param index1 The index of the parameter. 
	 * @return A pointer to the histogram 
	 */ 
	TH1D * MCMCGetH1Marginalized(int index1)
		{ return fMCMCH1Marginalized[index1]; }; 

	/*
	 * Returns a histogram of the 2-d marginalized distribution of the
	 * ith and jth parameter. 
	 * @param index1 The index of the first parameter. 
	 * @param index2 The index of the second parameter. 
	 * @return A pointer to the histogram. 
	 */ 
	TH2D * MCMCGetH2Marginalized(int index1, int index2); 

	/* @} */ 

	/** \name Setters */ 
	/* @{ */ 

	/*
	 * Sets the scale factor for the width of the trial function. 
	 * @param scale The scale factor. 
	 */ 
	void MCMCSetTrialFunctionScale(double scale) 
	{ fMCMCTrialFunctionScale = scale; }; 

	/*
	 * Sets the number of parameters of the Markov chain. 
	 * @param n The number of parameters. 
	 */ 
	void MCMCSetNParameters(int n); 

	/*
	 * Sets the number of Markov chains which are run in parallel.
	 * @param n The number of chains. 
	 */ 
	void MCMCSetNChains(int n); 

	/*
	 * Sets the maximum number of iterations. 
	 * @param n The maximum number of iterations. 
	 */ 
	void MCMCSetNIterationsMax(int n) 
	{ fMCMCNIterationsMax = n; }; 

	/*
	 * Sets the number of iterations. 
	 * @param n The number of iterations. 
	 */ 
	void MCMCSetNIterationsRun(int n) 
	{ fMCMCNIterationsRun = n; }; 

	/*
	 * Sets the number of iterations needed for burn-in. 
	 * @param n The number of iterations needed for burn-in. 
	 */ 
	void MCMCSetNIterationsBurnIn(int n)
	{ fMCMCNIterationsBurnIn = n; }; 

	/*
	 * Sets the number of iterations needed for PCA. 
	 * @param n The number of iterations needed for PCA. 
	 */ 
	void MCMCSetNIterationsPCA(int n)
	{ fMCMCNIterationsPCA = n; }; 
	
	/*
	 * Sets flag to automatically calculate the number of iterations of
	 * a Markov chain.
	 * @param flag The flag. 
	 */
	void MCMCSetIterationsAuto(bool flag) 
	{ fMCMCFlagIterationsAuto = flag; }; 

	/*
	 * Sets the minimum efficiency required for a chain
	 */ 
	void MCMCSetMinimumEfficiency(double efficiency) 
	{ fMCMCEfficiencyMin = efficiency; }; 

	/*
	 * Sets the maximum efficiency required for a chain
	 */ 
	void MCMCSetMaximumEfficiency(double efficiency) 
	{ fMCMCEfficiencyMax = efficiency; }; 

	/*
	 * Sets flag to write Markov chains to file.  
	 * @param flag The flag. 
	 */ 
	void MCMCSetWriteChainToFile(bool flag) 
	{ fMCMCFlagWriteChainToFile = flag; }; 

	/*
	 * Sets the initial position for a chain 
	 * @param chain The chain index. 
	 * @param x0 The intial position. 
	 * @see MCMCSetIntitialPositions. 
	 */ 
	void MCMCSetInitialPosition(std::vector<double> x0); 
	void MCMCSetInitialPosition(int chain, std::vector<double> x0); 

	/* 
	 * Sets the initial positions for all chains 
	 * @param x0s The initial positions for all chains. 
	 */ 
	void MCMCSetInitialPositions(std::vector<double> x0s); 
	
	/*
	 * Sets flag which defined initial position. 
	 */ 
	void MCMCSetFlagInitialPosition(int flag) 
	{ fMCMCFlagInitialPosition = flag; };  

	/*
	 * Sets the R-value criterion for convergence of all chains 
	 */ 
	void MCMCSetRValueCriterion(double r)
	{ fMCMCRValueCriterion = r; }; 

	/*
	 * Sets the flag to either perform a pre-run with PCA or not 
	 */ 
	void MCMCSetFlagPCA(bool flag)
	{ fMCMCFlagPCA = flag; }; 

	/*
	 * Sets the tree containing the Markov chains 
	 */ 
	void MCMCSetMarkovChainTrees(std::vector <TTree *> trees); 

	/*
	 * Set a flag to control if during the PCA the least eigenvectors
	 * should be ignored or not. 
	 */ 
	void MCMCSetFlagPCATruncate(bool flag) 
	{ fMCMCFlagPCATruncate = flag; }; 
	
	/*
	 * Sets the minimum ratio of an eigenvalue to the largest eigenvalue
	 * below which it is ignored if fMCMCFlagPCATruncate is true. 
	 */ 
	void MCMCSetPCAMinimumRatio(double ratio)
	{ fMCMCPCAMinimumRatio = ratio; }; 

	/* @} */ 

	/** \name Miscellaneous methods */ 
	/* @{ */ 

	/*
	 * Adds a parameter 
	 * @param min The minimum value of the parameter. 
	 * @param max The maximum value of the parameter. 
	 * @return The number of parameters. 
	 */ 
	int MCMCAddParameter(double min, double max); 

	/*
	 * Random walk trial function. The function is symmetric and
	 * used for the Metropolis algorithm. 
	 * @param x A point with the dimension
	 * fMCMCNParameters.
	 */ 
	void MCMCTrialFunction(std::vector <double> &x); 
	void MCMCTrialFunctionSingle(int ichain, int iparameter, std::vector <double> &x); 

	/*
	 * Independent chain trial function. The function does not
	 * have to be symmetric and is used for the
	 * Metropolis-Hastings algorithm.
	 * @param x A point with the dimension fMCMCNParameters. 
	 * @return The transition probability. 
	 */ 
	double MCMCTrialFunctionIndependent(std::vector <double> &xnew, std::vector <double> &xold, bool newpoint); 

	/*
	 * Trial function. 
	 * @param x A point with the dimension fMCMCNParameters. 
	 */ 
	void MCMCTrialFunctionAuto(std::vector <double> &x); 

	/*
	 * Trial function for the MCMC relative to the old point. No PCA. 
	 * @param pxold A pointer to the old point. The length of the vector 
	 * equals fMCMCNParameters.
	 * @param pxnew A pointer to the new point. The length of the vector 
	 * equals fMCMCNParameters. 
	 * @return flag_compute A flag which indicates whether to compute a 
	 * new point (true) or to just pass the value of the function (flase). 
	 * @return The value of the trial function. 
	 */ 
	double MCMCTrialFunctionRelativeNoPCA(std::vector <double> * xold, std::vector<double> * xnew, bool flag_compute); 

	/*
	 *
	 */
	void MCMCGetProposalPoint(int chain, std::vector <double> xnew, std::vector <double> xold); 

	/*
	 * Returns a trial point for the Metropolis algorithm. 
	 * @param chain The chain index. 
	 * @param x A proposal point. 
	 * @param pca A bool whether to use PCA or not 
	 * @return A flag indicating whether the new point lies within the allowed range. 
	 */ 
	bool MCMCGetProposalPointMetropolis(int chain, std::vector <double> &x, bool pca); 

	/*
	 * Returns a trial point for the Metropolis algorithm. 
	 * @param chain The chain index. 
	 * @param x A proposal point. 
	 * @param pca A bool whether to use PCA or not 
	 * @return A flag indicating whether the new point lies within the allowed range. 
	 */ 
	bool MCMCGetProposalPointMetropolis(int chain, int parameter, std::vector <double> &x, bool pca); 

	/*
	 * Returns a trial point for the Metropolis algorithm. 
	 * @param chain The chain index. 
	 * @param x A proposal point. 
	 * @param pca A bool whether to use PCA or not 
	 * @return A flag indicating whether the new point lies within the allowed range. 
	 */ 
	bool MCMCGetProposalPointMetropolisHastings(int chain, std::vector <double> &xnew, std::vector <double> &xold); 


	/* 
	 * Returns a new point using the Metropolis algorithm. 
	 * @param chain The chain index. 
	 * @param pca A bool whether to use PCA or not 
	 */ 
	bool MCMCGetNewPointMetropolis(int chain = 0, bool pca = false); 
	bool MCMCGetNewPointMetropolis(int chain = 0, int parameter = 0, bool pca = false); 

	/* 
	 * Returns a new point using the Metropolis algorithm. 
	 * @param chain The chain index. 
	 */ 
	bool MCMCGetNewPointMetropolisHastings(int chain = 0); 

	/* 
	 * Returns a new point using simulated annealing.
	 * @param chain The chain index. 
	 * @param pca A bool whether to use PCA or not 
	 */ 
	bool MCMCGetNewPointSimulatedAnnealing(int chain = 0, bool pca = false); 

	/*
	 * Returns the temperature accoring to the annealing schedule.
	 * @param chain The chain index. 
	 */
	double MCMCAnnealingSchedule(int chain); 

	/*
	 * Updates statistics, plots, etc. 
	 */ 
	void MCMCUpdateStatistics(); 

	void MCMCUpdateStatisticsCheckMaximum(); 

	void MCMCUpdateStatisticsFillHistograms();

	void MCMCUpdateStatisticsTestConvergenceAllChains(); 

	void MCMCUpdateStatisticsWriteChains(); 

	/*
	 * The probability density 
	 */ 
	virtual double LogEval(std::vector <double> parameters); 

	/*
	 * Perform a run for the PCA 
	 */ 
	void MCMCPCARun(); 

	/*
	 * Perform Metropolis algorithm 
	 */ 
	int MCMCMetropolis(); 

	/*
	 * Perform a pre run for the Metropolis algorithm 
	 */ 
	int MCMCMetropolisPreRun(); 

	/*
	 * Perform Metropolis-Hastings algorithm 
	 */ 
	int MCMCMetropolisHastings(); 

	/*
	 * Perform simulated annealing algorithm 
	 */ 
	int MCMCSimulatedAnnealing(); 

	/*
	 * Performs an initial run 
	 */ 
	void MCMCInitialRun(); 

	/*
	 * Resets the run statistics 
	 */ 
	void MCMCResetRunStatistics(); 

	/*
	 * Initializes the Markov chains
	 */ 
	void MCMCInitializeMarkovChains(); 

	/*
	 * Initializes the engine. 
	 */ 
	int MCMCInitialize(); 

	/*
	 * An interface to fitting functions. This function will be overloaded in BCModel. 
	 */
	virtual void MCMCUpdateFunctionFitting()
	  { return; }; 

	/*
	 * Allows user to interface during the run 
	 */ 
	virtual void MCMCUserInterface(); 

	/* @} */ 

 private: 

	/* 
	 * Copies this BCEngineMCMC into another one. 
	 */ 
	void Copy(BCEngineMCMC & enginemcmc) const; 

	/*
	 * Define a pointer to a member function 
	 */ 
	typedef bool (BCEngineMCMC::*MCMCPointerToGetProposalPoint) (int chain, std::vector <double> xnew, std::vector <double> xold) const; 

	/*
	 * Pointer to a member function 
	 */ 
	MCMCPointerToGetProposalPoint fMCMCPointerToGetProposalPoint; 

 protected: 

	/* 
	 * Number of parameters
	 */ 
	int fMCMCNParameters; 

	/* 
	 * Parameter boundaries 
	 */ 
	std::vector <double> fMCMCBoundaryMin; 
	std::vector <double> fMCMCBoundaryMax; 

	/*
	 * Number of Markov chains ran in parallel. 
	 */ 
	int fMCMCNChains;  

	/*
	 * Number of total iterations of the Markov chains. The length of
	 * the vector is equal to fMCMCNChains. 
	 */ 
	std::vector<int> fMCMCNIterations; 

	/*
	 * Number of iterations for updating scale factors 
	 */
	int fMCMCNIterationsUpdate; 

	/* 
	 * Number of iterations needed for each chain to convergence. The
	 * length of the vector is equal to fMCMCNChains. 
	 */ 
	std::vector<int> fMCMCNIterationsConvergenceLocal; 

	/* 
	 * Number of iterations needed for all chains to convergence simulaneously. 
	 */ 
	int fMCMCNIterationsConvergenceGlobal; 

	/*
	 * Flag for convergence 
	 */ 
	bool fMCMCFlagConvergenceGlobal;

	/*
	 * Maximum number of iterations for a Markov chain prerun. 
	 */  
	int fMCMCNIterationsMax; 

	/*
	 * Number of iterations for a Markov chain run. 
	 */  
	int fMCMCNIterationsRun; 

	/* 
	 * Number of iterations for burn-in. These iterations are not
	 * included in fMCMCNIterations. 
	 */
	int fMCMCNIterationsBurnIn; 

	/*
	 * Number of iterations for PCA. These iterations are not included
	 * in FMCMCNIterations.  
	 */ 
	int fMCMCNIterationsPCA; 

	/*
	 * Mean values of the data in the PCA coordinate system 
	 */ 
	std::vector <double> fMCMCPCAMean; 

	/*
	 * Variance of the data in the PCA coordinate system 
	 */ 
	std::vector <double> fMCMCPCAVariance; 

	/*
	 * Flag to control if during the PCA the least eigenvectors should
	 * be ignored or not.
	 */ 
	bool fMCMCFlagPCATruncate; 

	/*
	 * Minimum ratio of an eigenvalue to the largest eigenvalue below
	 * which it is ignored if fMCMCFlagPCATruncate is true.
	 */ 
	double fMCMCPCAMinimumRatio; 

	/*
	 * If the least eigenvectors are ignored this is the number of
	 * dimensions remaining.
	 */ 
	int fMCMCNDimensionsPCA; 

	/*
	 * Number of accepted trials and not accepted trials for each
	 * chain. The length of the vectors is equal to fMCMCNChains *
	 * fMCMCNParameters. For each chain these numbers add up to
	 * fMCMCNIterations.
	 */ 
	std::vector<int> fMCMCNTrialsTrue; 
	std::vector<int> fMCMCNTrialsFalse; 

	/* 
	 * The mean and variance of all values of each Markov chain. The
	 * length of the vectors is equal to fMCMCNChains.
	 */ 
	std::vector <double> fMCMCMean; 
	std::vector <double> fMCMCVariance; 

	/* 
	 * The sum and sum squared of all values of each Markov chain. These
	 * are used to calculate the mean and variance of each chain. The
	 * length of the vectors is equal to fMCMCNChains.
	 */ 
	std::vector <double> fMCMCSum; 
	std::vector <double> fMCMCSum2; 

	/*
	 * Flag to automatically calculate the number of iterations of a
	 * Markov chain. 
	 */ 
	bool fMCMCFlagIterationsAuto; 

	/*
	 * Flag to write Markov chains to file. 
	 */ 
	bool fMCMCFlagWriteChainToFile; 

	/*
	 * Scales the width of the trial functions by a global factor.
	 */ 
	double fMCMCTrialFunctionScale; 

	/*
	 * Scales the width of the trial functions by a scale factor
	 * for each parameter and chain
	 */ 
	std::vector <double> fMCMCTrialFunctionScaleFactor; 

	/* 
	 * Defines if a prerun has been performed or not 
	 */ 
	bool fMCMCFlagPreRun; 

	/*
	 * The intial position of each Markov chain. The length of the
	 * vectors is equal to fMCMCNChains * fMCMCNParameters. First, the
	 * values of the first Markov chain are saved, then those of the
	 * second and so on.
	 */
	std::vector <double> fMCMCInitialPosition; 

	/*
	 * Minimum efficiency for MCMC
	 */ 
	double fMCMCEfficiencyMin; 

	/*
	 * Maximum efficiency for MCMC
	 */ 
	double fMCMCEfficiencyMax; 

	/*
	 * Variable which defines the initial position. 0 (default) center
	 * of the allowed region, (1) random initial position (2)
	 * pre-defined intial position. 
	 */ 
	int fMCMCFlagInitialPosition; 

	/*
	 * The current points of each Markov chain. The length of the
	 * vectors is equal to fMCMCNChains * fMCMCNParameters. First, the
	 * values of the first Markov chain are saved, then those of the
	 * second and so on.
	 */
	std::vector <double> fMCMCx; 

	/*
	 * A temporary vector for a single Markov chain. 
	 */ 
	std::vector <double> fMCMCxLocal; 

	/* 
	 * The log of the probability of the current points of each Markov
	 * chain. The length of the vectors is fMCMCNChains.
	 */ 
	std::vector<double> fMCMCLogProbx; 

	/*
	 * The maximum points of each Markov chain. The length of the vector
	 * is fMCMCNChains * fMCMCNParameters. First, the values of the
	 * first Markov chain are saved, then those of the second and so on.
	 */ 
	std::vector <double> fMCMCMaximumPoints; 

	/*
	 * The maximum (log) probability of each Markov chain. The length of
	 * the vector is fMCMCNChains. 
	 */ 
	std::vector <double> fMCMCMaximumLogProb; 

	/*
	 * The R-value criterion for convergence. 
	 */ 
	double fMCMCRValueCriterion; 

	/*
	 * The R-value at which the chains did converge 
	 */ 
	double fMCMCRValue; 

	/* 
	 * The starting temperature for simulated annealing
	 */ 
	double fMCMCSimulatedAnnealingT0; 

	/*
	 * Random number generator 
	 */ 
	TRandom3 * fMCMCRandom; 

	/*
	 * PCA
	 */ 
	TPrincipal * fMCMCPCA; 

	/*
	 * Flag to perform PCA or not 
	 */ 
	bool fMCMCFlagPCA; 

	/*
	 * Number of bins per dimension for the marginalized distributions. 
	 */ 
	int fMCMCH1NBins; 
	int fMCMCH2NBins; 

	/*
	 * An array of marginalized distributions. 
	 */ 
	std::vector <TH1D *> fMCMCH1Marginalized; 
	std::vector <TH2D *> fMCMCH2Marginalized; 

	/*
	 * Control plots. 
	 */ 
	TH1D * fMCMCH1RValue;     // R-value 
	TH1D * fMCMCH1Efficiency; // Efficiency of the Markov chain 
  
	/*
	 * The trees containing the Markov chains. The length of the vector
	 * is fMCMCNChains. 
	 */ 
	std::vector<TTree *> fMCMCTrees; 

}; 

// --------------------------------------------------------- 

#endif 
