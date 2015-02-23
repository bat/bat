#ifndef __BCENGINEMCMC__H
#define __BCENGINEMCMC__H

/**
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
 * Copyright (C) 2007-2013, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

// ---------------------------------------------------------

// #include "BCVariable.h"
#include "BCParameter.h"
#include "BCParameterSet.h"
#include "BCObservable.h"
#include "BCVariableSet.h"
#include "BCLog.h"

#include <vector>
#include <limits>

#include <TMatrixDSym.h>
#include <TMatrixD.h>

// ROOT classes
class TF1;
class TH1;
class TH2;
class TTree;
class TFile;
class TRandom3;

class BCH1D;
class BCH2D;
class BCVariable;

// ---------------------------------------------------------

class BCEngineMCMC {
	
public:

	/** \name Enumerators  */
	/** @{ */

	/** An enumerator for the status of a test. */
	enum Precision {
		kLow,
		kQuick,
		kMedium,
		kHigh,
		kVeryHigh
	};

	/** An enumerator for the phase of the Markov chain. */
	enum MCMCPhase {
		kMCMCPreRunEfficiencyCheck       = -1,
		kMCMCPreRunConvergenceCheck      = -2,
		kMCMCPreRunFulfillMinimum        = -3,
		kMCMCUnsetPhase                  =  0,
		kMCMCMainRun                     = +1
	};

	/** An enumerator for markov-chain position initialization. */
	enum MCMCInitialPosition {
		kMCMCInitContinue         = -1, ///< continue markov chain from last value
		kMCMCInitCenter           =  0, ///< select centers of parameter ranges
		kMCMCInitRandomUniform    =  1, ///< randomly distribute uniformly over parameter ranges
		kMCMCInitUserDefined      =  2, ///< initialize to user-provided points
		kMCMCInitRandomPrior      =  3, ///< randomly distribute according to factorized priors
		kMCMCInitRandomGaussPrior =  4  ///< randomly distribute according to Gaussian assumption from factorized priors
	};

	/** @} */
	/** \name Structs */
	/** @{ */
	
	/** A struct for holding statistical information about samples. */
	struct MCMCStatistics {
		MCMCStatistics(const BCEngineMCMC::MCMCStatistics &other);
		MCMCStatistics(unsigned n_par=0, unsigned n_obs=0);
		unsigned n_samples;
		std::vector<double> mean;
		std::vector<double> variance;
		std::vector<std::vector<double> > covariance;
		std::vector<double> minimum;
		std::vector<double> maximum;
		double probability_mean;
		double probability_variance;
		double probability_mode;
		std::vector<double> mode;
		unsigned n_samples_efficiency;
		std::vector<double> efficiency;

		/** clear all members.
		 * @param clear_mode Flag for clearing information about mode*/
		void Clear(bool clear_mode=true, bool clear_effieciency=true);

		/** init all members
		 * @param n_par number of parameters
		 * @param n_obs number of observables. */
		void Init(unsigned n_par, unsigned n_obs);

		/** reset all members
		 * @param reset_mode flag for resetting information about mode.
		 * @param reset_efficiency flag for resetting information about efficiencies. */
		void Reset(bool reset_mode=true, bool reset_efficiency=true);
				
		/** reset efficiencies */
		void ResetEfficiencies();

		/** assignment operator. */
		MCMCStatistics & operator  = (const MCMCStatistics& rhs);
		/** addition assignment operator. */
		MCMCStatistics & operator += (const MCMCStatistics& rhs);

		/** update statistics given new values.
		 * @param par Current parameter values.
		 * @param obs Current user-define observables values. */
		void Update(double prob, const std::vector<double> & par, const std::vector<double> & obs);
	};

	/** @} */
	/** \name Constructors and destructors */
	/** @{ */

	/**
	 * Default constructor. */
	BCEngineMCMC(const char * name = "model");

	/**
	 * Default copy constructor. */
	BCEngineMCMC(const BCEngineMCMC & enginemcmc);

	/**
	 * Read in MCMC constructor.
	 * @param filename Path of file holding model.
	 * @param name Name of model (file should contain TTree's [name]_mcmc and [name]_parameters.\n
	 * if empty string is given, properly matching TTrees are searched for in the file.
	 * @param reuseObservables Flag for whether to load observables for file (true; default) or to let user respecify observables.*/
	BCEngineMCMC(std::string filename, std::string name, bool reuseObservables=true);

	/**
	 * Destructor. */
	virtual ~BCEngineMCMC();

	/** @} */
	/** \name Assignment operators */
	/** @{ */

	/**
	 * Defaut assignment operator */
	BCEngineMCMC & operator = (const BCEngineMCMC & engineMCMC)
	{ Copy(engineMCMC); return *this; }
	
	/** @} */
	/** \name Getters */
	/** @{ */

	/**
	 * @return The name of the engine. */
	const std::string & GetName() const
	{ return fName; }

	/**
	 * @return The name of the engine with spaces removed. */
	const std::string & GetSafeName() const
	{ return fSafeName; }

	/**
	 * @return number of Markov chains */
	unsigned MCMCGetNChains() const
	{ return fMCMCNChains; }

	/**
	 * @return lag of the Markov chains */
	unsigned MCMCGetNLag() const
	{ return fMCMCNLag; }

	/**
	 * @return number of iterations */
	const std::vector<unsigned> & MCMCGetNIterations() const
	{ return fMCMCNIterations; }

	/**
	 * @return current iterations */
	int MCMCGetCurrentIteration() const
	{ return fMCMCCurrentIteration; }

	/**
	 * @return current chain index */
	int MCMCGetCurrentChain() const
	{ return fMCMCCurrentChain; }

	/**
	 * @return number of iterations needed for all chains to
	 * converge simultaneously */
	unsigned MCMCGetNIterationsConvergenceGlobal() const
	{ return fMCMCNIterationsConvergenceGlobal; }
	
	/**
	 * @return minimum number of pre-run iterations for a Markov chain */
	unsigned MCMCGetNIterationsPreRunMin() const
	{ return fMCMCNIterationsPreRunMin; }

	/**
	 * @return maximum number of pre-run iterations for a Markov chain */
	unsigned MCMCGetNIterationsPreRunMax() const
	{ return fMCMCNIterationsPreRunMax; }

	/**
	 * @return number of iterations for a Markov chain */
	unsigned MCMCGetNIterationsRun() const
	{ return fMCMCNIterationsRun; }

	/**
	 * @return number of iterations for an efficiency check. */
	unsigned MCMCGetNIterationsEfficiencyCheck() const
	{ return fMCMCNIterationsEfficiencyCheck; }

	/**
	 * @return number of iterations after statistics update. */
	unsigned MCMCGetNIterationsConvergenceCheck() const
	{ return fMCMCNIterationsConvergenceCheck; }

	/**
	 * @return number of iterations after statistics clear. */
	unsigned MCMCGetNIterationsClearConvergenceStats() const
	{ return fMCMCNIterationsClearConvergenceStats; }

	/**
	 * @return minimum efficiency required for a chain. */
	double MCMCGetMinimumEfficiency() const
	{ return fMCMCEfficiencyMin; }

	/**
	 * @return maximum efficiency required for a chain. */
	double MCMCGetMaximumEfficiency() const
	{ return fMCMCEfficiencyMax; }

	/**
	 * @return trial function scale factor lower limit */
	double MCMCGetTrialFunctionScaleFactorLowerLimit() const
	{ return fMCMCScaleFactorLowerLimit; }

	/**
	 * @return trial function scale factor upper limit */
	double MCMCGetTrialFunctionScaleFactorUpperLimit() const
	{ return fMCMCScaleFactorUpperLimit; }
    
	/**
	 * @return scale factor for all parameters and chains */
	const std::vector<std::vector<double> > & MCMCGetTrialFunctionScaleFactor() const
	{ return fMCMCTrialFunctionScaleFactor; }

	/**
	 * @return scale factor for all parameters and achain.
	 * @param c chain index */
	std::vector<double> MCMCGetTrialFunctionScaleFactor(unsigned c) const
	{ return (c<fMCMCTrialFunctionScaleFactor.size()) ? fMCMCTrialFunctionScaleFactor[c] : std::vector<double>(); }
	
	/**
	 * @return scale factor for a parameter and a chain.
	 * @param c chain index
	 * @param p parameter index */
	double MCMCGetTrialFunctionScaleFactor(unsigned c, unsigned p) const
	{ return (c<fMCMCTrialFunctionScaleFactor.size() and p<fMCMCTrialFunctionScaleFactor[c].size()) ? fMCMCTrialFunctionScaleFactor[c][p] : std::numeric_limits<double>::quiet_NaN(); }
	
	/**
	 * @return current point of each Markov chain */
	const std::vector<std::vector<double> > & MCMCGetx() const
	{ return fMCMCx; }

	/**
	 * @param c index of the Markov chain
	 * @return current point of the Markov chain */
	std::vector<double> MCMCGetx(unsigned c) const
	{ return (c<fMCMCx.size()) ? fMCMCx[c] : std::vector<double>(); }
	
	/**
	 * @param c chain index
	 * @param p parameter index
	 * @return parameter of the Markov chain */
	double MCMCGetx(unsigned c, unsigned p) const
	{ return (c<fMCMCx.size() and p<fMCMCx[c].size()) ? fMCMCx[c][p] : std::numeric_limits<double>::quiet_NaN(); }

	/**
	 * @return log of the probability of the current points of each Markov chain */
	const std::vector<double> & MCMCGetLogProbx() const
	{ return fMCMCprob; }
	
	/**
	 * @return log of the probability of the current points of the Markov chain.
	 * @param c chain index */
	double MCMCGetLogProbx(unsigned c) const
	{ return (c<fMCMCprob.size()) ? fMCMCprob[c] : -std::numeric_limits<double>::infinity(); }

	/**
	 * @return pointer to the phase of a run. */
	BCEngineMCMC::MCMCPhase MCMCGetPhase() const
	{ return fMCMCPhase; }

	/**
	 * @return flag which defined initial position */
	BCEngineMCMC::MCMCInitialPosition MCMCGetFlagInitialPosition() const
	{ return fMCMCFlagInitialPosition; }

	/**
	 * @return expansion factor to be used with initial-position option kMCMCInitRandomGaussPrior. */
	double MCMCGetInitialPositionExpansionFactor() const
	{ return fMCMCInitialPositionExpansionFactor; }

	/**
	 * @return whether to use a multivariate proposal function. */
	bool MCMCGetMultivariateProposalFunction() const
	{ return fMCMCMultivariateProposalFunction; }

	// /**
	//  * @return multivariate-proposal-function tuning schedule parameter. */
	// double MCMCGetMultivariateProposalFunctionTuningScheduleParameter() const
	// { return fMultivariateProposalFunctionTuningScheduleParameter; }

	/**
	 * @return multivariate-proposal-function Cholesky decomposition nudge size. */
	double MCMCGetMultivariateProposalFunctionEpsilon() const
	{ return fMultivariateProposalFunctionEpsilon; }

	/**
	 * @return multivariate-proposal-function scale multiplier. */
	double MCMCGetMultivariateProposalFunctionScaleMultiplier() const
	{ return fMultivariateProposalFunctionScaleMultiplier; }

	/**
	 * @return R-value criterion */
	double MCMCGetRValueCriterion() const
	{ return fMCMCRValueCriterion; }

	/**
	 * @return R-value criterion for parameters */
	double MCMCGetRValueParametersCriterion() const
	{ return fMCMCRValueParametersCriterion; }

	/**
	 * @return R-value */
	double MCMCGetRValue() const
	{ return fMCMCRValue; }

	/**
	 * @return R-value for a parameter
	 * @param i parameter index */
	double MCMCGetRValueParameters(unsigned i) const
	{ return (i<fMCMCRValueParameters.size()) ? fMCMCRValueParameters[i] : std::numeric_limits<double>::infinity(); }

	/** Flag for correcting convergence checking for initial sampling variability. */
	bool MCMCGetCorrectRValueForSamplingVariability() const
	{ return fCorrectRValueForSamplingVariability; }

	/**
	 * @return the flag if MCMC has been performed or not */
	bool MCMCGetFlagRun() const
	{ return fMCMCFlagRun; }

	/**
	 * @return the flag if MCMC pre-run should be performed or not */
	bool MCMCGetFlagPreRun() const
	{ return fMCMCFlagPreRun; }

	/**
	 * Retrieve the tree containing the Markov chain. */
	TTree * MCMCGetMarkovChainTree()
	{ return fMCMCTree;}

	/**
	 * Retrieve the tree containing the parameter information. */
	TTree * MCMCGetParameterTree()
	{ return fParameterTree;}

	/**
	 * Retrieve output file for MCMC. */
	TFile * MCMCGetOutputFile()
	{ return fMCMCOutputFile; }

	/**
	 * Get combined statistics for all chains. */
	BCEngineMCMC::MCMCStatistics MCMCGetStatistics() const
	{ return fMCMCStatistics_AllChains; }
	    
	/**
	 * Get MCMC statistics for one chain.
	 * @param c Chain to get statisics of. */
	BCEngineMCMC::MCMCStatistics MCMCGetStatistics(unsigned c) const
	{ return (c<fMCMCStatistics.size()) ? fMCMCStatistics[c] : BCEngineMCMC::MCMCStatistics(); }
	
	/**
	 * Close MCMC output file. */
	void MCMCCloseOutputFile();

	/**
	 * @return Flag for whether to rescale histogram ranges to fit MCMC reach after pre-run. */
	bool GetRescaleHistogramRangesAfterPreRun() const
	{ return fRescaleHistogramRangesAfterPreRun; }

	/**
	 * @return Factor by which to enlarge histogram ranges when rescaling to add padding beyond range. */
	double GetHistogramRescalePadding() const
	{ return fHistogramRescalePadding; }

	/**
	 * @param index Index of histogram of which to check existence
	 * @return Whether the marginalized histogram exists. */
	bool MarginalizedHistogramExists(unsigned index) const
	{ return index<fH1Marginalized.size() and fH1Marginalized[index]; }

	/**
	 * @param index1 X Index of histogram of which to check existence.
	 * @param index2 Y Index of histogram of which to check existence.
	 * @return Whether the marginalized histogram exists. */
	bool MarginalizedHistogramExists(unsigned index1, unsigned index2) const
	{ return index1<fH2Marginalized.size() and index2<fH2Marginalized[index1].size() and fH2Marginalized[index1][index2]; }

	/**
	 * Obtain the individual marginalized distributions
	 * with respect to one parameter as a ROOT TH1
	 * @note The most efficient method is to access by index.
	 * @param name The parameter's name
	 * @return 1D marginalized probability */
	TH1 * GetMarginalizedHistogram(const char * name) const
	{ return GetMarginalizedHistogram(fParameters.Index(name)); }
   
	/**
	 * Obtain the individual marginalized distributions
	 * with respect to one parameter as a ROOT TH1
	 * @param index The parameter index
	 * @return 1D marginalized probability */
	TH1 * GetMarginalizedHistogram(unsigned index) const;
   
	/**
	 * Obtain the individual marginalized distributions
	 * with respect to two parameters as a ROOT TH2.
	 * @note The most efficient method is to access by indices.
	 * @param name1 Name of first parameter
	 * @param name2 Name of second parameter
	 * @return 2D marginalized probability */
	TH2 * GetMarginalizedHistogram(const char * name1, const char * name2) const
	{ return GetMarginalizedHistogram(fParameters.Index(name1), fParameters.Index(name2)); }
   
	/**
	 * Obtain the individual marginalized distributions
	 * with respect to two parameters as a ROOT TH2.
	 * @param index1 Index of first parameter
	 * @param index2 Index of second parameter
	 * @return 2D marginalized probability */
	TH2 * GetMarginalizedHistogram(unsigned index1, unsigned index2) const;

	/**
	 * Obtain the individual marginalized distributions
	 * with respect to one parameter.
	 * @note The most efficient method is to access by index.
	 * @note Ownership of the returned heap object is conferred to the caller.
	 * @param name The parameter's name
	 * @return 1D marginalized probability */
	BCH1D * GetMarginalized(const char * name) const
	{ return GetMarginalized(fParameters.Index(name)); }
   
	/**
	 * Obtain the individual marginalized distributions
	 * with respect to one parameter.
	 * @note Ownership of the returned heap object is conferred to the caller.
	 * @param index The parameter index
	 * @return 1D marginalized probability */
	BCH1D * GetMarginalized(unsigned index) const;
   
	/**
	 * Obtain the individual marginalized distributions
	 * with respect to two parameters.
	 * @note The most efficient method is to access by indices.
	 * @note Ownership of the returned heap object is conferred to the caller.
	 * @param name1 Name of first parameter
	 * @param name2 Name of second parameter
	 * @return 2D marginalized probability */
	BCH2D * GetMarginalized(const char * name1, const char * name2) const
	{ return GetMarginalized(fParameters.Index(name1), fParameters.Index(name2)); }
   
	/**
	 * Obtain the individual marginalized distributions
	 * with respect to two parameters.
	 * @note Ownership of the returned heap object is conferred to the caller.
	 * @param index1 Index of first parameter
	 * @param index2 Index of second parameter
	 * @return 2D marginalized probability */
	BCH2D * GetMarginalized(unsigned index1, unsigned index2) const;
	
	/**
	 * @param observables Whether to check max length of user-defined observable names.
	 * @return length of longest parameter name. */
	unsigned GetMaximumParameterNameLength(bool observables=true) const
	{ return (observables) ? std::max(fParameters.MaxNameLength(),fObservables.MaxNameLength()) : fParameters.MaxNameLength(); }

	/**
	 * @param index The index of the observable running first over 
	 * 0,...,N_parameters in the ParameterSet, and then over
	 * N_parameters,...,N_parameters+N_observables in the ObservableSet
	 * @return The observable. */
	BCVariable * GetVariable(unsigned index) const
	{ return ((index<GetNParameters()) ? fParameters.Get(index) : ((index<GetNVariables()) ? fObservables.Get(index-GetNParameters()) : NULL)); }

	/**
	 * @return The number of parameters of the model. */
	unsigned GetNVariables() const
	{ return fParameters.Size() + fObservables.Size(); }

	/**
	 * @return Parameter set. */
	BCParameterSet & GetParameters()
	{ return fParameters; }

	/**
	 * @param index The index of the parameter in the parameter set.
	 * @return The parameter. */
	BCParameter * GetParameter(int index) const
	{ return dynamic_cast<BCParameter*>(fParameters.Get(index)); }

	/**
	 * @param name The name of the parameter in the parameter set.
	 * @return The parameter. */
	BCParameter * GetParameter(const char * name) const
	{ return dynamic_cast<BCParameter*>(fParameters.Get(name)); }

	/**
	 * @return The number of parameters of the model. */
	unsigned GetNParameters() const
	{ return fParameters.Size(); }

	/**
	 * @return The number of fixed parameters. */
	unsigned GetNFixedParameters() const
	{ return fParameters.GetNFixedParameters(); }

	/**
	 * @return The number of free parameters. */
	unsigned GetNFreeParameters() const
	{ return fParameters.GetNFreeParameters(); }

	/**
	 * @return Observable set. */
	BCVariableSet & GetObservables()
	{ return fObservables; }

	/**
	 * @param index The index of the observable in the observable set.
	 * @return The user-defined observable. */
	BCObservable * GetObservable(int index) const
	{ return dynamic_cast<BCObservable*>(fObservables.Get(index)); }

	/**
	 * @param name The name of the observable in the observable set.
	 * @return The user-defined observable. */
	BCObservable * GetObservable(const char * name) const
	{ return dynamic_cast<BCObservable*>(fObservables.Get(name)); }

	/**
	 * @return The number of user-defined observables. */
	unsigned GetNObservables() const
	{ return fObservables.Size(); }

	/**
	 * @return vector of parameter and observable values at global mode. */
	virtual const std::vector<double> & GetGlobalMode() const
	{ return fMCMCStatistics_AllChains.mode; }

	/**
	 * @return vector of parameter values at global mode. */
	virtual std::vector<double> GetBestFitParameters() const;

	/**
	 * @return vector of observable values at global mode. */
	virtual std::vector<double> GetBestFitObservables() const;

	/**
	 * @return vector of the local modes of parameters and observables
	 * @param force_recalculation flag for forcing recalculation of local modes from histograms. */
	const std::vector<double> & GetLocalModes(bool force_recalculation=false);

	/**
	 * @return The log of the value at the mode. */
	virtual double GetLogMaximum() const
	{ return fMCMCStatistics_AllChains.probability_mode; };

	/** 
	 * @return Flag whether to reuse user-defined observables from MCMC tree when looping through it. */ 
	bool MCMCGetReuseObservables() const
	{ return fMCMCTreeReuseObservables; }

	/**
	 * @return BCH1D object that stores drawing options for all BCH1D's. */
	BCH1D * GetBCH1DdrawingOptions()
	{ return fBCH1DdrawingOptions; }

	/**
	 * @return BCH2D object that stores drawing options for all BCH2D's. */
	BCH2D * GetBCH2DdrawingOptions()
	{ return fBCH2DdrawingOptions; }

	/** @} */
	/** \name Setters */
	/** @{ */

	/**
	 * Sets the name of the engine.
	 * @param name Name of the engine */
	void SetName(const char * name);

	/**
	 * Sets the name of the engine.
	 * @param name Name of the engine */
	void SetName(const std::string name)
	{ SetName(name.data()); }

	/**
	 * Set scale factor lower limit */
	void MCMCSetScaleFactorLowerLimit(double l)
	{ fMCMCScaleFactorLowerLimit = l; }

	/**
	 * Set scale factor upper limit */
	void MCMCSetScaleFactorUpperLimit(double l)
	{ fMCMCScaleFactorUpperLimit = l; }

	/**
	 * Set the scale factors for the trial functions
	 * @param scale a vector of doubles containing the scale factors */
	void MCMCSetTrialFunctionScaleFactor(std::vector<double> scale)
	{ fMCMCTrialFunctionScaleFactorStart = scale; }

	/**
	 * Set flag to automatically set initial scale factors by factorized prior variances. */
	void MCMCAutoSetTrialFunctionScaleFactors(bool flag=true)
	{ fMCMCAutoSetTrialFunctionScaleFactors = flag; }

	/**
	 * Sets the number of Markov chains which are run in parallel. */
	void MCMCSetNChains(unsigned n)
	{ fMCMCNChains = n; }

	/**
	 * Sets the lag of the Markov chains */
	void MCMCSetNLag(unsigned n)
	{ fMCMCNLag = n; }

	/**
	 * Sets current chain number. */
	void MCMCSetCurrentChain(int n)
	{ fMCMCCurrentChain = n; }

	/**
	 * Sets the maximum number of iterations in the pre-run. */
	void MCMCSetNIterationsPreRunMax(unsigned n)
	{ fMCMCNIterationsPreRunMax = n; }

	/**
	 * Sets the number of iterations. */
	void MCMCSetNIterationsRun(unsigned n)
	{ fMCMCNIterationsRun = n; }

	/**
	 * Sets the minimum number of iterations in the pre-run */
	void MCMCSetNIterationsPreRunMin(unsigned n)
	{ fMCMCNIterationsPreRunMin = n; }

	/**
	 * Sets the number of iterations in the pre-run after which the
	 * efficiency is checked and proposal scales are updated. */
	void MCMCSetNIterationsEfficiencyCheck(unsigned n)
	{ fMCMCNIterationsEfficiencyCheck = n; }

	/**
	 * Sets the number of iterations in the pre-run after which an
	 * check on the convergence is done.
	 * @param n The number of iterations.*/
	void MCMCSetNIterationsConvergenceCheck(unsigned n)
	{ fMCMCNIterationsConvergenceCheck = n; }

	/**
	 * Sets the number of iterations in the pre-run after which the
	 * convergence data is cleared.
	 * @param n The number of iterations.*/
	void MCMCSetNIterationsClearConvergenceStats(unsigned n)
	{ fMCMCNIterationsClearConvergenceStats = n; }

	/**
	 * Sets the minimum efficiency required for a chain. */
	void MCMCSetMinimumEfficiency(double efficiency)
	{ fMCMCEfficiencyMin = efficiency; }

	/**
	 * Sets the maximum efficiency required for a chain. */
	void MCMCSetMaximumEfficiency(double efficiency)
	{ fMCMCEfficiencyMax = efficiency; }

	/**
	 * Set the random number seed */
	void MCMCSetRandomSeed(unsigned seed);

	/**
	 * Sets flag to write Markov chains to file. */
	void MCMCSetWriteChainToFile(bool flag)
	{ fMCMCFlagWriteChainToFile = flag; }

	/**
	 * Sets flag to write pre run to file. */
	void MCMCSetWritePreRunToFile(bool flag)
	{ fMCMCFlagWritePreRunToFile = flag; }

	/**
	 * Sets the initial positions for all chains.
	 * @param x0s initial positions for all chains. */
	void MCMCSetInitialPositions(const std::vector<double> &x0s);

	/**
	 * Sets the initial positions for all chains.
	 * @param x0s initial positions for all chains. */
	void MCMCSetInitialPositions(const std::vector< std::vector<double> > & x0s)
	{ fMCMCInitialPosition = x0s; MCMCSetFlagInitialPosition(BCEngineMCMC::kMCMCInitUserDefined); }

	/**
	 * Sets flag which defines initial position.  */
	void MCMCSetFlagInitialPosition(BCEngineMCMC::MCMCInitialPosition flag)
	{ fMCMCFlagInitialPosition = flag; }

	/**
	 * Sets expansion factor for use with initial position setting kMCMCInitRandomGaussPrior. */
	void MCMCSetInitialPositionExpansionFactor(double f)
	{ fMCMCInitialPositionExpansionFactor = fabs(f); }

	/**
	 * Sets the flag which controls the sequence parameters during the
	 * running of the MCMC.  */
	void MCMCSetMultivariateProposalFunction(bool flag)
	{ fMCMCMultivariateProposalFunction = flag; }

	// /**
	//  * Sets multivariate-proposal-function covariance update tuning parameter. */
	// void MCMCSetMultivariateProposalFunctionTuningScheduleParameter(double lambda)
	// { fMultivariateProposalFunctionTuningScheduleParameter = lambda; }

	/**
	 * Sets multivariate-proposal-function cholesky-decomposition nudge. */
	void MCMCSetMultivariateProposalFunctionEpsilon(double epsilon)
	{ fMultivariateProposalFunctionEpsilon = epsilon; }

	/**
	 * Sets multivariate-proposal-function scale multiplier. */
	void MCMCSetMultivariateProposalFunctionScaleMultiplier(double s)
	{ fMultivariateProposalFunctionScaleMultiplier = s; }

	/** Sets whether to fill histograms.
	 * Applies to 1D and 2D histogams of all parameters and observables so far added. */
	void MCMCSetFlagFillHistograms(bool flag)
	{ fParameters.FillHistograms(flag); fObservables.FillHistograms(flag); }

	/** Sets the whether to fill histograms.
	 * Applies to all parameters and observables so far added. */
	void MCMCSetFlagFillHistograms(bool flag_1d, bool flag_2d)
	{ fParameters.FillHistograms(flag_1d,flag_2d); fObservables.FillHistograms(flag_1d,flag_2d); }

	/** Sets the flag if a prerun should be performed or not. */
	void MCMCSetFlagPreRun(bool flag)
	{ fMCMCFlagPreRun = flag; }

	/**
	 * Sets the R-value criterion for convergence of all chains. */
	void MCMCSetRValueCriterion(double r)
	{ fMCMCRValueCriterion = r; }

	/**
	 * Sets the parameter R-value criterion for convergence of all chains */
	void MCMCSetRValueParametersCriterion(double r)
	{ fMCMCRValueParametersCriterion = r; }

	/** Set flag to correct convergence checking for initial sampling variability. */
	void MCMCSetCorrectRValueForSamplingVariability(bool flag=true)
	{ fCorrectRValueForSamplingVariability = flag; }

	/**
	 * Initialize the trees containing the Markov chains and parameter info. */
	virtual void InitializeMarkovChainTree(bool replacetree=false, bool replacefile=false);

	/**
	 * Set the precision for the MCMC run. */
	void MCMCSetPrecision(BCEngineMCMC::Precision precision);

	/**
	 * Copy precision for the MCMC run from other model. */
	void MCMCSetPrecision(const BCEngineMCMC * other);
	
	/**
	 * Set the number of bins for the marginalized distribution of all parameters
	 * @param nbins Number of bins */
	void SetNbins(unsigned int nbins)
	{ fParameters.SetNBins(nbins); fObservables.SetNBins(nbins); }

	/** 
	 * @param flag Flag whether to reuse user-defined observables from MCMC tree when looping through it. */ 
	void MCMCSetReuseObservables(bool flag)
	{ fMCMCTreeReuseObservables = flag; }

	/**
	 * Set flag for rescaling histogram ranges after pre-run. */
	void SetRescaleHistogramRangesAfterPreRun(bool flag)
	{ fRescaleHistogramRangesAfterPreRun = flag; }

	/**
	 * Set enlargement factor of range for when rescaling. */
	void SetHistogramRescalingPadding(double factor)
	{ fHistogramRescalePadding = factor; }

	/**
	 * Turn on/off writing of Markov chain to root file.
	 * If setting true, use function with filename arguments.
	 * @param flag Flag for writing Markov chain to ROOT file (true) or not (false). */
	void WriteMarkovChain(bool flag);

	/** Turn on writing of Markov chain to root file.
	 * @param filename Name of file to write chain to.
	 * @param option file-open options (TFile), must be "NEW", "CREATE", "RECREATE", or "UPDATE" (i.e. writeable). */
	void WriteMarkovChain(std::string filename, std::string option);

	/**
	 * Write marginalization histograms to file. */
	void WriteMarginalizedDistributions(std::string filename, std::string option);

	/** @} */
	/** \name Prior setting functions (all deprecated).
	 * \brief The priors are not used for the likelihood calculation by
	 * BCEngineMCMC, but are used for initializing the positions of the
	 * chains. */
	/** @{ */

	/**
	 * @deprecated Instead call: GetParameter(index)->SetPriorConstant()
	 * Set constant prior for this parameter
	 * @param index the index of the parameter
	 * @return success of action. */
	bool SetPriorConstant(unsigned index);

	/**
	 * @deprecated Instead call: GetParameter(name)->SetPriorConstant()
	 * Set constant prior for this parameter
	 * @param name the name of the parameter
	 * @return Success of action */
	bool SetPriorConstant(const char* name)
	{ return SetPriorConstant(fParameters.Index(name)); }

	/**
	 * @deprecated Instead call: GetPrarameter(index)->SetPrior(new BCTF1Prior(f))
	 * Set prior for a parameter.
	 * @param index The parameter index
	 * @param f A pointer to a function describing the prior
	 * @return success of action. */
	bool SetPrior(unsigned index, TF1* const f);

	/**
	 * @deprecated Instead call: GetParameter(name)->SetPrior(new BCTF1Prior(f))
	 * Set prior for a parameter.
	 * @param name The parameter name
	 * @param f A pointer to a function describing the prior
	 * @return success of action. */
	bool SetPrior(const char* name, TF1* const f)
	{ return SetPrior(fParameters.Index(name), f); }

	/**
	 * @deprecated Instead call: GetParameter(index)->Fix(value)
	 * Fixes parameter to value.
	 * @param index The parameter index
	 * @param value The position of the delta function.
	 * @return success of action. */
	bool SetPriorDelta(unsigned index, double value);

	/**
	 * @deprecated Instead call: GetParameter(name)->Fix(value)
	 * Fixes parameter to value.
	 * @param name The parameter name
	 * @param value The position of the delta function.
	 * @return success of action. */
	bool SetPriorDelta(const char* name, double value)
	{ return SetPriorDelta(fParameters.Index(name),value); }

	/**
	 * @deprecated Instead call: GetParameter(index)->SetPrior(new BCGaussianPrior(mean,sigma))
	 * Set Gaussian prior for a parameter.
	 * @param index The parameter index
	 * @param mean The mean of the Gaussian
	 * @param sigma The sigma of the Gaussian
	 * @return success of action. */
	bool SetPriorGauss(unsigned index, double mean, double sigma);

	/**
	 * @deprecated Instead call: GetParameter(name)->SetPrior(new BCGaussianPrior(mean,sigma))
	 * Set Gaussian prior for a parameter.
	 * @param name The parameter name
	 * @param mean The mean of the Gaussian
	 * @param sigma The sigma of the Gaussian
	 * @return success of action. */
	bool SetPriorGauss(const char* name, double mean, double sigma)
	{ return SetPriorGauss(fParameters.Index(name), mean, sigma); }

	/**
	 * @deprecated Instead call: GetParameter(index)->SetPrior(new BCSplitGaussianPrior(mean,sigma_below,sigma_above))
	 * Set Gaussian prior for a parameter with two different widths.
	 * @param index The parameter index
	 * @param mean The mean of the Gaussian
	 * @param sigma_below Standard deviation below mean.
	 * @param sigma_above Standard deviation above mean.
	 * @return success of action. */
	bool SetPriorGauss(unsigned index, double mean, double sigma_below, double sigma_above);

	/**
	 * @deprecated Instead call: GetParameter(name)->SetPrior(new BCSplitGaussianPrior(mean,sigma_below,sigma_above))
	 * Set Gaussian prior for a parameter with two different widths.
	 * @param name The parameter name
	 * @param mean The mean of the Gaussian
	 * @param sigmadown The sigma (down) of the Gaussian
	 * @param sigmaup The sigma (up)of the Gaussian
	 * @return success of action. */
	int SetPriorGauss(const char* name, double mean, double sigmadown, double sigmaup)
	{	return SetPriorGauss(fParameters.Index(name), mean, sigmadown, sigmaup); }

	/**
	 * @deprecated Instead call: GetParameter(index)->SetPrior(new BCTH1Prior(h,interpolate))
	 * Set prior for a parameter.
	 * @param index parameter index
	 * @param h pointer to a histogram describing the prior
	 * @param interpolate whether or not to use linear interpolation
	 * @return success of action. */
	bool SetPrior(unsigned index, TH1 * const h, bool interpolate=false);

	/**
	 * @deprecated Instead call: GetParameter(name)->SetPrior(new BCTH1Prior(h,interpolate))
	 * Set prior for a parameter.
	 * @param name parameter name
	 * @param h pointer to a histogram describing the prior
	 * @param interpolate whether or not to use linear interpolation
	 * @return success of action. */
	bool SetPrior(const char* name, TH1 * const h, bool interpolate=false)
	{ return SetPrior(fParameters.Index(name),h,interpolate); }

	/**
	 * @deprecated Instead call: GetParameters().SetPriorConstantAll()
	 * Enable caching the constant value of the prior, so LogAPrioriProbability
	 * is called only once. Note that the prior for ALL parameters is
	 * assumed to be constant. The value is computed from
	 * the parameter ranges, so make sure these are defined before this method is
	 * called. */
	void SetPriorConstantAll()
	{ fParameters.SetPriorConstantAll(); }

	/** @} **/
	/** \name Miscellaneous methods */
	/** @{ */

	/**
	 * Create histograms from parameter and observable sets.
	 * @param rescale_ranges Rescale axis ranges to range reached by MCMC tree if true */
	virtual void CreateHistograms(bool rescale_ranges=false);

	/**
	 * Prints a summary on the screen. */
	virtual void PrintSummary();

	/**
	 * Prints a summary to a file.
	 * @return Success of action.*/
	virtual bool PrintResults(const char * file) const;

	/**
	 * Print parameters
	 * @P is a vector of the parameter values to be printed
	 * @output is a pointer to the output function to be used,
	 * which defaults to BCLog::OutSummary */
	void PrintParameters(std::vector<double> const & P, void (*output)(const char *) = BCLog::OutSummary) const;

	/**
	 * Print all marginalizations.
	 * @param filename Path to file to print to
	 * @param hdiv Number of columns of plots per page
	 * @param vdiv Number of rows of plots per page
	 * @return Number of plots printed */
	unsigned PrintAllMarginalized(std::string filename, unsigned hdiv=1, unsigned vdiv=1) const;
	
	/**
	 * Print a summary plot for the parameters and user-defined observables.
	 * @param filename Path to filename to print to.
	 * @param npar Number of parameters per page, print all on one page if set to zero or negative
	 * @param interval_content Probability mass to display in smallest X interval band
	 * @param quantile_values Vector of quantile values to draw
	 * @param rescale_ranges Flag for rescaling to range surveyed by MCMC chains
	 * @return Number of pages printed. */
	unsigned PrintParameterPlot(std::string filename, int npar=10, double interval_content=68e-2, std::vector<double> quantile_vals=std::vector<double>(0), bool rescale_ranges=true) const;
			
	/**
	 * Draw a summary plot for the parameters in the range provided to current pad
	 * @par i0 Index of first parameter to print.
	 * @par npar Number of parameters to print, set to 0 to print all.
	 * @param interval_content Probability mass to display in smallest X interval band
	 * @param quantile_values Vector of quantile values to draw
	 * @param rescale_ranges Flag for rescaling to range surveyed by MCMC chains
	 * @return Success of action. */
	bool DrawParameterPlot(unsigned i0, unsigned npar=0, double interval_content=68e-2, std::vector<double> quantile_vals=std::vector<double>(0), bool rescale_ranges=true) const;
			
	/**
	 * Print a correlation matrix for the parameters.
	 * @return An error flag. */
	int PrintCorrelationMatrix(const char * filename = "matrix.pdf");
			
	/**
	 * Print a correlation plot for the parameters.
	 * @return An error flag. */
	int PrintCorrelationPlot(const char * filename = "correlation.pdf", bool include_observables=true);
			
	/**
	 * Print a Latex table of the parameters.
	 * @return An error flag. */
	int PrintParameterLatex(const char * filename);

	/**
	 * Copy object
	 * @param enginemcmc Object to copy from */
	void Copy(const BCEngineMCMC & enginemcmc);

	/**
	 * Adds a parameter.
	 * @param min minimum value of the parameter
	 * @param max maximum value of the parameter
	 * @param latexname Optional latexname used for plotting
	 * @return number of parameters after adding */
	virtual int AddParameter(const char * name, double min, double max, const char * latexname = "")
	{ return AddParameter(new BCParameter(name, min, max, latexname)); }

	/**
	 * Adds a parameter to the model.
	 * @param parameter A model parameter
	 * @see AddParameter(const char * name, double lowerlimit, double upperlimit, const char * latexname); */
	virtual int AddParameter(BCParameter* parameter)
	{ return fParameters.Add(parameter); }

	/**
	 * Adds a user-calculated observable.
	 * @param min minimum value of the observable
	 * @param max maximum value of the observable
	 * @param latexname Optional latexname used for plotting
	 * @return number of observables after adding */
	virtual int AddObservable(const char * name, double min, double max, const char * latexname = "")
	{ return AddObservable(new BCObservable(name,min,max,latexname)); }

	/**
	 * Adds a user-calculated observable to the model.
	 * @param observable A user-calculated observable
	 * @see AddObservable(const char * name, double lowerlimit, double upperlimit, ObservableFunction * fn, const char * latexname); */
	virtual int AddObservable(BCObservable* obs)
	{ return fObservables.Add(obs); }

	/**
	 * Evaluates user-defined observables at current state of all
	 * chains and stores results in fMCMCObservables*/
	virtual void EvaluateObservables();

	/**
	 * Evaluates user-defined observables at current state of chain
	 * and stores results in fMCMCObservables*/
	virtual void EvaluateObservables(unsigned chain);

	/**
	 * Evaluates user-defined observables To be overloaded by user
	 * to calculate user-observables. */
	virtual void CalculateObservables(const std::vector<double> &/*pars*/)
	{}

	/**
	 * Random walk trial function. The default trial function is a
	 * Breit-Wigner. It can be overloaded by the user to set the trial
	 * function.
	 * @param ichain the chain index
	 * @param x point with the dimension fMCMCNParameters */
	virtual void MCMCTrialFunction(unsigned ichain, std::vector<double> &x);

	/**
	 * Random walk trial function. The default trial function is a
	 * Breit-Wigner. It can be overloaded by the user to set the trial
	 * function.
	 * @param ichain the chain index
	 * @param ipar the parameter index
	 * @return the unscaled proposal point */
	virtual double MCMCTrialFunctionSingle(unsigned ichain, unsigned ipar);

	/**
	 * Returns a trial point for the Metropolis algorithm.
	 * @param chain chain index
	 * @param x proposal point
	 * @return flag indicating whether the new point lies within the allowed range */
	bool MCMCGetProposalPointMetropolis(unsigned chain, std::vector<double> &x);

	/**
	 * Returns a trial point for the Metropolis algorithm.
	 * @param chain chain index
	 * @param x proposal point
	 * @return flag indicating whether the new point lies within the allowed range */
	bool MCMCGetProposalPointMetropolis(unsigned chain, unsigned parameter, std::vector<double> &x);

	/**
	 * Generates a new point using the Metropolis algorithm.
	 * @param chain chain index */
	bool MCMCGetNewPointMetropolis();
	bool MCMCGetNewPointMetropolis(unsigned chain);
	bool MCMCGetNewPointMetropolis(unsigned chain, unsigned parameter);

	/**
	 * Updates statistics: fill marginalized distributions */
	void MCMCInChainFillHistograms();

	// /**
	//  * Updates statistics: check convergence */
	// void MCMCInChainTestConvergenceAllChains();

	/**
	 * Updates statistics: write chains to file */
	void MCMCInChainWriteChains();

	/**
	 * Needs to be overloaded in the derived class.
	 * @return natural logarithm of the function to map with MCMC */
	virtual double LogEval(const std::vector<double> & parameters) = 0;

	/**
	 * Runs Metropolis algorithm.
	 * @return Success of action. */
	bool MCMCMetropolis();

	/**
	 * Runs a pre run for the Metropolis algorithm.
	 * @return Success of action. */
	bool MCMCMetropolisPreRun();

	/**
	 * Initializes the engine.
	 * @return Success of action. */
	bool MCMCInitialize();

	/**
	 * Reset the MCMC variables.
	 */
	virtual void ResetResults();

	/**
	 * Empty the sequence of parameters.
	 *
	 * @param hard Delete the parameters. @warning Could lead to problems if multiple modes use identical parameters.
	 */
	void ClearParameters(bool hard=false)
	{ fParameters.Clear(hard); }

	/**
	 * Interface allowing to execute arbitrary code for each iteration
	 * of the MCMC. The frequency of calling this method is influenced
	 * by the setup of the Lag and whether or not the MCMC is run with
	 * ordered parameters. This method needs to be overloaded in the derived
	 * class. */
	virtual void MCMCIterationInterface() 
	{ MCMCUserIterationInterface(); }

	/**
	 * Method executed for every iteration of the MCMC. User's code should be
	 * provided via overloading in the derived class*/
	virtual void MCMCUserIterationInterface()
	{}

	/**
	 * Interface allowing to execute arbitrary code for each new point
	 * of the MCMC. This method needs to be overloaded in the derived
	 * class
	 * @param point point that was generated and checked
	 * @param ichain index of the chain
	 * @param accepted flag whether or not the point was accepted for the chain
	 */
	virtual void MCMCCurrentPointInterface(std::vector<double> & /*point*/, int /*ichain*/, bool /*accepted*/)
	{}

	/**
	 * Load parameters and observables from tree. */
	virtual bool LoadParametersFromTree(TTree * partree, bool reuseObservables=true);

	/**
	 * Check parameter tree against model. */
	virtual bool ParameterTreeMatchesModel(TTree * partree, bool checkObservables=true);

	/**
	 * Load previous MCMC run.
	 * @param filename Pathname of file containing model trees.
	 * @param mcmcTreeName Name of tree inside file containing MCMC, empty string (default) loads [modelname]_mcmc.
	 * @param parameterTreeName Name of tree inside file containing parameter list, empty string (default) loads [modelname]_parameters.
	 * @param reuseObservables Flag for whether to load observables from parameter list and MCMC trees. */
	virtual bool LoadMCMC(std::string filename, std::string mcmcTreeName="", std::string parameterTreeName="", bool reuseObservables=true);

	/**
	 * Load previous MCMC run. */
	virtual bool LoadMCMC(TTree * mcmcTree, TTree * parTree, bool reuseObservables=true);

	/**
	 * Check tree structure for MCMC tree. */
	virtual bool ValidMCMCTree(TTree * tree, bool checkObservables=true) const;

	/**
	 * Check tree structure for parameter tree. */
	virtual bool ValidParameterTree(TTree * tree) const;

	/**
	 * Marginalize from TTree. */
	virtual void Remarginalize(bool autorange=true);

	/**
	 * Calculcate evidnce. */
	virtual double CalculateEvidence(double epsilon=1e-2);

	/**
	 * Update cholesky-decompositions for multivariate proposal functio. */
	virtual bool UpdateCholeskyDecompositions();

	/** @} */

private:


	/**
	 * Defines a type of a pointer to a member function. */
	typedef bool (BCEngineMCMC::*MCMCPointerToGetProposalPoint) (int chain, std::vector<double> xnew, std::vector<double> xold) const;

	/**
	 * Pointer to a member function */
	MCMCPointerToGetProposalPoint fMCMCPointerToGetProposalPoint;

	/**
	 * Keeps variables that need to be both members and thread local.
	 */
	struct MCMCThreadLocalStorage
	{
		/**
		 * Store local proposal point
		 */
		std::vector<double> xLocal;
				
		/**
		 * Random number generator
		 */
		TRandom3 * rng;

		/**
		 * ctor
		 * @param dim Dimension of a temporary parameter vector
		 */
		MCMCThreadLocalStorage(const unsigned & dim);

		MCMCThreadLocalStorage(const MCMCThreadLocalStorage & other);

		MCMCThreadLocalStorage & operator = (const MCMCThreadLocalStorage & other);

		virtual ~MCMCThreadLocalStorage();
	};

	/**
	 * Keep thread local variables private.
	 */
	std::vector<MCMCThreadLocalStorage> fMCMCThreadLocalStorage;

	/**
	 * Ensure that there are as many storages as chains
	 */
	void SyncThreadStorage();

protected:

	/**
	 * Print model summary to stream. */
	virtual void PrintSummaryToStream(std::ofstream & ofi) const;

	/**
	 * Print best fit to stream. */
	virtual void PrintBestFitToStream(std::ofstream & ofi) const;

	/**
	 * Print marginalization to stream. */
	virtual void PrintMarginalizationToStream(std::ofstream & ofi) const;

	/**
	 * Update Paramater TTree with scales and efficiencies. */
	void UpdateParameterTree();

	/**
	 * Name of the engine. */
	std::string fName;

	/**
	 * Safe name of the engine for use in naming ROOT objects. */
	std::string fSafeName;

	/**
	 * return appropriate update interval */
	unsigned  UpdateFrequency(unsigned N);
    
	/**
	 * Parameter settings */
	BCParameterSet fParameters;

	/**
	 * User-calculated Observables Set */
	BCVariableSet fObservables;

	/**
	 * Number of Markov chains ran in parallel */
	unsigned fMCMCNChains;

	/**
	 * The lag for the Markov Chain */
	unsigned fMCMCNLag;

	/**
	 * Number of total iterations of the Markov chains. The length of
	 * the vector is equal to fMCMCNChains. */
	std::vector<unsigned> fMCMCNIterations;

	/**
	 * The current iteration number. If not called within the running
	 * of the algorithm, return -1. */
	int fMCMCCurrentIteration;

	/**
	 * The current chain index. If not called within the running of the
	 * algorithm, return -1. */
	int fMCMCCurrentChain;

	/**
	 * Number of iterations for checking efficiencies and updating scale
	 * factors */
	unsigned fMCMCNIterationsEfficiencyCheck;

	/**
	 * Number of iterations for updating scale factors */
	unsigned fMCMCNIterationsConvergenceCheck;

	/**
	 * Number of iterations for clearing data for updating scale factors */
	unsigned fMCMCNIterationsClearConvergenceStats;

	/**
	 * Number of iterations needed for all chains to convergence
	 * simultaneously */
	int fMCMCNIterationsConvergenceGlobal;

	/**
	 * Maximum number of iterations for a Markov chain prerun */
	unsigned fMCMCNIterationsPreRunMax;

	/**
	 * Number of iterations for a Markov chain run */
	unsigned fMCMCNIterationsRun;

	/**
	 * Minimum number of iterations for the pre-run */
	unsigned fMCMCNIterationsPreRunMin;

	/**
	 * Flag to write Markov chains to file */
	bool fMCMCFlagWriteChainToFile;

	/**
	 * Flag to write pre run to file */
	bool fMCMCFlagWritePreRunToFile;

	/*
	 * Output file for for writing MCMC Tree. */
	TFile * fMCMCOutputFile;

	/*
	 * Output filename for for writing MCMC Tree. */
	std::string fMCMCOutputFilename;

	/*
	 * Output file open option for for writing MCMC Tree. */
	std::string fMCMCOutputFileOption;

	/**
	 * Lower limit for scale factors */
	double fMCMCScaleFactorLowerLimit;
    
	/**
	 * Upper limit for scale factors */
	double fMCMCScaleFactorUpperLimit;

	/**
	 * Scale factors for proposal functions.
	 * First index is for the chain,
	 * second index is for the parameters (if not using Multivariate proposal function). */
	std::vector<std::vector<double> > fMCMCTrialFunctionScaleFactor;

	/**
	 * Start values of the scale factors for the trial functions. */
	std::vector<double> fMCMCTrialFunctionScaleFactorStart;

	/**
	 * flag to auto set initial trial function scale factors. */
	bool fMCMCAutoSetTrialFunctionScaleFactors;

	/**
	 * Covariance matrices used in multivariate proposal functions. */
	std::vector<TMatrixDSym> fMultivariateProposalFunctionCovariance;

	/**
	 * Cholesky decompositions for multivariate proposal function.
	 * Index is over chain. */
	std::vector<TMatrixD> fMultivariateProposalFunctionCholeskyDecomposition;

	/**
	 * number of multivariate-proposal-function tuning steps performed. */
	unsigned fMultivariateProposalFunctionTuningSteps;

	// /**
	//  * multivariate-proposal-function tuning-schedule parameter. */
	// double fMultivariateProposalFunctionTuningScheduleParameter;

	/**
	 * multivariate-proposal-function cholesky-decomposition nudge. */
	double fMultivariateProposalFunctionEpsilon;

	/**
	 * factor to multiply or divide scale factors by in adjusting multivariate-proposal-function scales. */
	double fMultivariateProposalFunctionScaleMultiplier;

	/**
	 * Defines if a prerun should be performed or not */
	bool fMCMCFlagPreRun;

	/**
	 * Defines if MCMC has been performed or not */
	bool fMCMCFlagRun;

	/**
	 * The intial position of each Markov chain. The length of the
	 * vectors is equal to fMCMCNChains * fMCMCNParameters. First, the
	 * values of the first Markov chain are saved, then those of the
	 * second and so on */
	std::vector<std::vector<double> > fMCMCInitialPosition;

	/**
	 * Factor to multiply Gaussian approximations of prior standard deviations
	 * when generating initial points accordingly. */
	double fMCMCInitialPositionExpansionFactor;

	/**
	 * The efficiencies for all parameters and chains. */
	std::vector<std::vector<double> > fMCMCEfficiencies;

	/**
	 * The minimum required efficiency for MCMC */
	double fMCMCEfficiencyMin;

	/**
	 * The maximum allowed efficiency for MCMC */
	double fMCMCEfficiencyMax;

	/**
	 * Variable which defines the initial position.
	 * See enum MCMCInitialPosition for possible values. */
	BCEngineMCMC::MCMCInitialPosition fMCMCFlagInitialPosition;

	/**
	 * Flag for using multivariate proposal function. */
	bool fMCMCMultivariateProposalFunction;

	/**i
	 * The phase of the run.
	 */
	BCEngineMCMC::MCMCPhase fMCMCPhase;

	/**
	 * The current points of each Markov chain. */
	std::vector<std::vector<double> > fMCMCx;
	
	/**
	 * The current values of the user-defined observables for each
	 * Markov chain. */
	std::vector<std::vector<double> > fMCMCObservables;

	/** 
	 * Statistics for each Markov chain. */
	std::vector<BCEngineMCMC::MCMCStatistics> fMCMCStatistics;

	/**
	 * Statistics across all Markov chains. */
	BCEngineMCMC::MCMCStatistics fMCMCStatistics_AllChains;

	/**
	 * The log of the probability of the current points of each Markov
	 * chain. The length of the vectors is fMCMCNChains. */
	std::vector<double> fMCMCprob;

	/**
	 * The log of the likelihood of the current points of each Markov
	 * chain. The length of the vectors is fMCMCNChains. */
	std::vector<double> fMCMCLogLikelihood;

	/**
	 * The log of the likelihood of the proposed points of each Markov
	 * chain. The length of the vectors is fMCMCNChains. */
	std::vector<double> fMCMCLogLikelihood_Provisional;

	/**
	 * The log of the prior of the current points of each Markov
	 * chain. The length of the vectors is fMCMCNChains. */
	std::vector<double> fMCMCLogPrior;

	/**
	 * The log of the prior of the proposed points of each Markov
	 * chain. The length of the vectors is fMCMCNChains. */
	std::vector<double> fMCMCLogPrior_Provisional;

	/**
	 * flag for correcting R value for initial sampling variability.
	 */
	bool fCorrectRValueForSamplingVariability;

	/**
	 * The R-value criterion for convergence of log-likelihood*/
	double fMCMCRValueCriterion;

	/**
	 * The R-value criterion for convergence of parameters */
	double fMCMCRValueParametersCriterion;

	/**
	 * The R-value at which the chains did converge */
	double fMCMCRValue;

	/** The R-values for each parameter */
	std::vector<double> fMCMCRValueParameters;

	/**
	 * Random number generator */
	TRandom3 * fRandom;

	/**
	 * An array of marginalized distributions */
	std::vector<TH1*>               fH1Marginalized;
	std::vector<std::vector<TH2*> > fH2Marginalized;

	/**
	 * The tree containing the Markov chains.*/
	TTree * fMCMCTree;

	/**
	 * flag for whether MCMC Tree successfully loaded.*/
	bool fMCMCTreeLoaded;
	
	/**
	 * flag for whether to reuse MCMC Tree's observables. */
	bool fMCMCTreeReuseObservables;

	/**
	 * fMCMCTree's variables. */
	unsigned int fMCMCTree_Chain;
	unsigned int fMCMCTree_Iteration;
	double fMCMCTree_Prob;
	double fMCMCTree_LogLikelihood;
	double fMCMCTree_LogPrior;
	std::vector<double> fMCMCTree_Parameters;
	std::vector<double> fMCMCTree_Observables;

	/**
	 * The tree containing the parameter information.*/
	TTree * fParameterTree;
	
	/**
	 * Vector of local modes. */
	std::vector<double> fLocalModes;

	/**
	 * A BCH1D (with no histogram) for storing BCH1D drawing options. */
	BCH1D * fBCH1DdrawingOptions;
	    
	/**
	 * A BCH2D (with no histogram) for storing BCH2D drawing options. */
	BCH2D * fBCH2DdrawingOptions;

	/**
	 * flag for rescaling of histograms after pre-run. */
	bool fRescaleHistogramRangesAfterPreRun;

	/**
	 * factor for enlarging range of histograms when rescaling. */
	double fHistogramRescalePadding;

};

// ---------------------------------------------------------

#endif
