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
  BCEngineMCMC(int n = 0); 

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
   * Returns the current point of each Markov chain. 
   */ 
  std::vector <double> MCMCGetx()
    { return fMCMCx; }; 

  /* 
   * Returns the current point of the ith Markov chain.
   */ 
  std::vector <double> * MCMCGetP2x() 
		{ return &fMCMCx; }; 

  /* 
   * Returns the current point of the ith Markov chain.
   * @param i The index of the Markov chain. 
   */ 
  std::vector <double> MCMCGetx(int i); 

  /* 
   * Returns the log of the probability of the current points of each
   * Markov chain. 
   */ 
  std::vector <double> MCMCGetLogProbx()
    { return fMCMCLogProbx; }; 

	/* 
   * Returns the pointer to the log of the probability of the current
   * points of each Markov chain.
   */ 
  std::vector <double> * MCMCGetP2LogProbx()
    { return &fMCMCLogProbx; }; 

  /*
   * Returns the minimum points of each Markov chain. 
   */ 
  std::vector <double> MCMCGetMinimumPoints()
    { return fMCMCMinimumPoints; }; 

  /*
   * Returns the minimum point of each ith Markov chain. 
   * @param i The index of the Markov chain. 
   */ 
  std::vector <double> MCMCGetMinimumPoints(int i); 

  /*
   * Returns the minimum (log) probability of each Markov chain.
   */ 
  std::vector <double> MCMCGetMinimumLogProb()
    { return fMCMCMinimumLogProb; }; 

  /*
   * Returns the trees containing the Markov chains. 
   */ 
  std::vector <TTree *> MCMCGetTrees()
    { return fMCMCTrees; }; 

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
	 * Trial function. 
	 * @param x A point with the dimension fMCMCNParameters. 
	 */ 
	void MCMCTrialFunction(std::vector <double> &x); 

	/*
	 * Trial function. 
	 * @param x A point with the dimension fMCMCNParameters. 
	 */ 
	void MCMCTrialFunctionAuto(std::vector <double> &x); 

	/*
	 * Returns a trial point. 
	 * @param chain The chain index. 
	 * @param x A proposal point. 
	 * @param pca A bool whether to use PCA or not 
	 * @return A flag indicating whether the new point lies within the allowed range. 
	 */ 
	bool MCMCGetProposalPoint(int chain, std::vector <double> &x, bool pca); 
	
	/* 
	 * Returns a new point using the Metropolis algorithm. 
	 * @param chain The chain index. 
	 * @param pca A bool whether to use PCA or not 
	 */ 
	bool MCMCGetNewPointMetropolis(int chain = 0, bool pca = false); 

	/*
	 * Updates the convergence criteria 
	 */ 
	void MCMCUpdateConvergence(); 

	/*
	 * The probability density 
	 */ 
	double LogEval(std::vector <double> parameters); 

	/*
	 * Perform a run for the PCA 
	 */ 
	void MCMCPCARun(); 

	/*
	 * Perform Metropolis algorithm 
	 */ 
	int MCMCMetropolis(); 

	/*
	 * Performs an initial run 
	 */ 
	void MCMCInitialRun(); 

	/*
	 * Resets the run statistics 
	 */ 
	void MCMCResetRunStatistics(); 

  /*
   * Initializes the engine. 
   */ 
  int MCMCInitialize(); 

	/*
	 * Allows user to interface during the run 
	 */ 
	void MCMCUserInterface(); 

	/* @} */ 

private: 

  /* 
   * Copies this BCEngineMCMC into another one. 
   */ 
  void Copy(BCEngineMCMC & enginemcmc) const; 

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
   * Number of iterations needed for each chain to convergence. The
   * length of the vector is equal to fMCMCNChains. 
   */ 
  std::vector<int> fMCMCNIterationsConvergenceLocal; 

  /* 
   * Number of iterations needed for all chains to convergence simulaneously. 
   */ 
  int fMCMCNIterationsConvergenceGlobal; 

  /*
   * Maximum number of iterations for a Markov chain. 
   */  
  int fMCMCNIterationsMax; 

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
   * Number of accepted trials and not accepted trials for each
   * chain. The length of the vectors is equal to fMCMCNChains. For
   * each chain these numbers add up to fMCMCNIterations. 
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
   * The intial position of each Markov chain. The length of the
   * vectors is equal to fMCMCNChains * fMCMCNParameters. First, the
   * values of the first Markov chain are saved, then those of the
   * second and so on.
   */
  std::vector <double> fMCMCInitialPosition; 

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
   * The minimum points of each Markov chain. The length of the vector
   * is fMCMCNChains * fMCMCNParameters. First, the values of the
   * first Markov chain are saved, then those of the second and so on.
   */ 
  std::vector <double> fMCMCMinimumPoints; 

  /*
   * The minimum (log) probability of each Markov chain. The length of
   * the vector is fMCMCNChains. 
   */ 
  std::vector <double> fMCMCMinimumLogProb; 

	/*
	 * The R-value criterion for convergence. 
	 */ 
  double fMCMCRValueCriterion; 

	/*
	 * The R-value at which the chains did converge 
	 */ 
  double fMCMCRValue; 

  /*
   * The trees containing the Markov chains. The length of the vector
   * is fMCMCNChains. 
   */ 
  std::vector<TTree *> fMCMCTrees; 

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

}; 

// --------------------------------------------------------- 

#endif 
