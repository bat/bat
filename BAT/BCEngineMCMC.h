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

#include "BCVariable.h"
#include "BCParameter.h"
#include "BCParameterSet.h"
#include "BCObservable.h"
#include "BCObservableSet.h"
#include "BCLog.h"

#include <vector>

class BCH1D;
class BCH2D;

// ROOT classes
class TH1D;
class TH2D;
class TTree;
class TRandom3;

// ---------------------------------------------------------

class BCEngineMCMC
{

   public:

      /** \name Enumerators  */
      /** @{ */

      /** An enumerator for the status of a test. */
	    enum Precision{ kLow, kQuick, kMedium, kHigh, kVeryHigh };

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
       * Destructor. */
      virtual ~BCEngineMCMC();

      /** @} */
      /** \name Assignment operators */
      /** @{ */

      /**
       * Defaut assignment operator */
      BCEngineMCMC & operator = (const BCEngineMCMC & engineMCMC);

      /** @} */
      /** \name Getters */
      /** @{ */

      /**
       * @return The name of the engine. */
      const std::string & GetName() const
         { return fName; }

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
      unsigned MCMCGetCurrentIteration() const
         { return fMCMCCurrentIteration; }

      /**
       * @return current chain index */
      unsigned MCMCGetCurrentChain() const
         { return fMCMCCurrentChain; }

      /**
       * @return number of iterations needed for all chains to
       * converge simultaneously */
      unsigned MCMCGetNIterationsConvergenceGlobal() const
         { return fMCMCNIterationsConvergenceGlobal; }

      /**
       * @return maximum number of iterations for a Markov chain */
      unsigned MCMCGetNIterationsMax() const
         { return fMCMCNIterationsMax; }

      /**
       * @return number of iterations for a Markov chain */
      unsigned MCMCGetNIterationsRun() const
         { return fMCMCNIterationsRun; }

      /**
       * @return minimum number of pre-run iterations for a Markov chain */
      unsigned MCMCGetNIterationsPreRunMin() const
         { return fMCMCNIterationsPreRunMin; }

      /**
       * @return number of iterations for an efficiency check. */
      unsigned MCMCGetNIterationsEfficiencyCheck() const
         { return fMCMCNIterationsEfficiencyCheck; }

      /**
       * @return number of iterations after statistics update. */
      unsigned MCMCGetNIterationsUpdate() const
         { return fMCMCNIterationsUpdate; }

      /**
       * @return number of iterations after statistics clear. */
      unsigned MCMCGetNIterationsUpdateClear() const
         { return fMCMCNIterationsUpdateClear; }

      /**
       * @return minimum efficiency required for a chain. */
      double MCMCGetMinimumEfficiency() const
	       { return fMCMCEfficiencyMin; }

      /**
       * @return maximum efficiency required for a chain. */
      double MCMCGetMaximumEfficiency() const
         { return fMCMCEfficiencyMax; }

      /**
       * @returns number of accepted trials for each parameter of each chain */
	    std::vector<std::vector<int> > MCMCGetNTrialsTrue() const
	       { return fMCMCNTrialsTrue; }

      /**
       * @returns number of trials */
      int MCMCGetNTrials() const
         { return fMCMCNTrials; }

      /**
       * @return mean value of the probability for each chain up to
       * the current iteration  */
      const std::vector<double> & MCMCGetprobMean() const
         { return fMCMCprobMean; }

      /**
       * @return mean value of the probability for each chain up to
       * the current iteration */
      const std::vector<double> & MCMCGetVariance() const
         { return fMCMCprobVar; }

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
       * @param ichain chain index */
      std::vector<double> MCMCGetTrialFunctionScaleFactor(unsigned ichain) const;

      /**
       * @return scale factor for a parameter and a chain.
       * @param ichain chain index
       * @param ipar parameter index */
      double MCMCGetTrialFunctionScaleFactor(unsigned ichain, unsigned ipar);

      /**
       * @return current point of each Markov chain */
	    const std::vector<std::vector<double> > & MCMCGetx() const
         { return fMCMCx; }

      /**
       * @param ichain index of the Markov chain
       * @return current point of the Markov chain */
      std::vector<double> MCMCGetx(unsigned ichain);

      /**
       * @param ichain chain index
       * @param ipar parameter index
       * @return parameter of the Markov chain */
      double MCMCGetx(unsigned ichain, unsigned ipar) const;

      /**
       * @return log of the probability of the current points of each Markov chain */
      const std::vector<double> & MCMCGetLogProbx() const
         { return fMCMCprob; }

      /**
       * @return log of the probability of the current points of the Markov chain.
       * @param ichain chain index */
      double MCMCGetLogProbx(unsigned ichain);

      /**
       * @return pointer to the phase of a run. */
      int MCMCGetPhase() const
         { return fMCMCPhase; }

      /**
       * @return maximum points of each Markov chain */
    	const std::vector<std::vector<double> > & MCMCGetMaximumPoints() const
         { return fMCMCxMax; }

      /**
       * @return maximum point of  Markov chain
       * @param i The index of the Markov chain */
      std::vector<double> MCMCGetMaximumPoint(unsigned i) const;

      /**
       * @return maximum (log) probability of each Markov chain */
      const std::vector<double> & MCMCGetMaximumLogProb() const
         { return fMCMCprobMax; }

      /**
       * @return flag which defined initial position */
      int MCMCGetFlagInitialPosition() const
         { return fMCMCFlagInitialPosition; }

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
      double MCMCGetRValueParameters(unsigned i)
         { return fMCMCRValueParameters.at(i); }

      /** Use strict or relaxed rule for Gelman/Rubin R-value */
      bool MCMCGetRValueStrict() const
         { return fMCMCRValueUseStrict; }

      /**
       * @return the flag if MCMC has been performed or not */
      bool MCMCGetFlagRun() const
         { return fMCMCFlagRun; }

      /**
       * Retrieve the tree containing the Markov chain.
       * @param i index of the Markov chain
       * @return pointer to the tree */
      TTree * MCMCGetMarkovChainTree(unsigned i)
      { if (i < fMCMCTrees.size())
          return fMCMCTrees.at(i);
        else
          return 0; }

	    /**
			 * @return The size of vector of marginalized 1D histograms */
	    unsigned GetN1DMarginalizations()
	       { return fH1Marginalized.size(); }
		 
	    /**
			 * @return The size of vector of vector of marginalized 2D histograms */
	    unsigned GetN2DMarginalizations()
	       { return fH2Marginalized.size(); }

	    /**
			 * @param i index of vector of which the size should be returned
			 * @return The size of the i'th vector of marginalized 2D histograms */
	    unsigned GetN2DMarginalizations(unsigned i)
	       { return (i<GetN2DMarginalizations()) ? fH2Marginalized[i].size() : 0; }

	    /**
			 * @param index Index of histogram of which to check existence
			 * @return Whether the marginalized histogram exists. */
	    bool MarginalizedHistogramExists(unsigned index)
	       { return index<fH1Marginalized.size() and fH1Marginalized[index]; }

	    /**
			 * @param index1 X Index of histogram of which to check existence.
			 * @param index2 Y Index of histogram of which to check existence.
			 * @return Whether the marginalized histogram exists. */
	    bool MarginalizedHistogramExists(unsigned index1, unsigned index2)
	       { return index1<fH2Marginalized.size() and index2<fH2Marginalized[index1].size() and fH2Marginalized[index1][index2]; }

      /**
       * Obtain the individual marginalized distributions
       * with respect to one parameter as a ROOT TH1D
       * @note The most efficient method is to access by index.
       * @param parameter Model parameter
       * @return 1D marginalized probability */
	    TH1D * GetMarginalizedHistogram(const BCParameter * parameter) const
	       { return GetMarginalizedHistogram(fParameters.Index(parameter->GetName())); }
   
      /**
       * Obtain the individual marginalized distributions
       * with respect to one parameter as a ROOT TH1D
       * @note The most efficient method is to access by index.
       * @param name The parameter's name
       * @return 1D marginalized probability */
      TH1D * GetMarginalizedHistogram(const char * name) const
         { return GetMarginalizedHistogram(fParameters.Index(name)); }
   
      /**
       * Obtain the individual marginalized distributions
       * with respect to one parameter as a ROOT TH1D
       * @param index The parameter index
       * @return 1D marginalized probability */
      TH1D * GetMarginalizedHistogram(unsigned index) const;
   
      /**
       * Obtain the individual marginalized distributions
       * with respect to two parameters as a ROOT TH2D.
       * @note The most efficient method is to access by indices.
       * @param parameter1 First parameter
       * @param parameter2 Second parameter
       * @return 2D marginalized probability */
      TH2D * GetMarginalizedHistogram(const BCParameter * parameter1, const BCParameter * parameter2) const
	       { return GetMarginalizedHistogram(fParameters.Index(parameter1->GetName()),fParameters.Index(parameter2->GetName())); }
   
      /**
       * Obtain the individual marginalized distributions
       * with respect to two parameters as a ROOT TH2D.
       * @note The most efficient method is to access by indices.
       * @param name1 Name of first parameter
       * @param name2 Name of second parameter
       * @return 2D marginalized probability */
      TH2D * GetMarginalizedHistogram(const char * name1, const char * name2) const
         { return GetMarginalizedHistogram(fParameters.Index(name1), fParameters.Index(name2)); }
   
      /**
       * Obtain the individual marginalized distributions
       * with respect to two parameters as a ROOT TH2D.
       * @param index1 Index of first parameter
       * @param index2 Index of second parameter
       * @return 2D marginalized probability */
      TH2D * GetMarginalizedHistogram(unsigned index1, unsigned index2) const;

      /**
       * Obtain the individual marginalized distributions
       * with respect to one parameter.
       * @note The most efficient method is to access by index.
       * @note Ownership of the returned heap object is conferred to the caller.
       * @param parameter Model parameter
       * @return 1D marginalized probability */
	    BCH1D * GetMarginalized(const BCParameter * parameter)
	       { return GetMarginalized(fParameters.Index(parameter->GetName())); }
   
      /**
       * Obtain the individual marginalized distributions
       * with respect to one parameter.
       * @note The most efficient method is to access by index.
       * @note Ownership of the returned heap object is conferred to the caller.
       * @param name The parameter's name
       * @return 1D marginalized probability */
      BCH1D * GetMarginalized(const char * name)
         { return GetMarginalized(fParameters.Index(name)); }
   
      /**
       * Obtain the individual marginalized distributions
       * with respect to one parameter.
       * @note Ownership of the returned heap object is conferred to the caller.
       * @param index The parameter index
       * @return 1D marginalized probability */
      BCH1D * GetMarginalized(unsigned index);
   
      /**
       * Obtain the individual marginalized distributions
       * with respect to two parameters.
       * @note The most efficient method is to access by indices.
       * @note Ownership of the returned heap object is conferred to the caller.
       * @param parameter1 First parameter
       * @param parameter2 Second parameter
       * @return 2D marginalized probability */
      BCH2D * GetMarginalized(const BCParameter * parameter1, const BCParameter * parameter2)
	       { return GetMarginalized(fParameters.Index(parameter1->GetName()),fParameters.Index(parameter2->GetName())); }
   
      /**
       * Obtain the individual marginalized distributions
       * with respect to two parameters.
       * @note The most efficient method is to access by indices.
       * @note Ownership of the returned heap object is conferred to the caller.
       * @param name1 Name of first parameter
       * @param name2 Name of second parameter
       * @return 2D marginalized probability */
      BCH2D * GetMarginalized(const char * name1, const char * name2)
         { return GetMarginalized(fParameters.Index(name1), fParameters.Index(name2)); }
   
      /**
       * Obtain the individual marginalized distributions
       * with respect to two parameters.
       * @note Ownership of the returned heap object is conferred to the caller.
       * @param index1 Index of first parameter
       * @param index2 Index of second parameter
       * @return 2D marginalized probability */
      BCH2D * GetMarginalized(unsigned index1, unsigned index2);
   
	    /**
			 * @param observables Whether to check max length of user-defined observable names.
			 * @return length of longest parameter name. */
 	    unsigned GetMaximumParameterNameLength(bool observables=true)
	       { return (observables) ? std::max(fParameters.MaxNameLength(),fObservables.MaxNameLength()) : fParameters.MaxNameLength(); }

      /**
       * @param index The index of the observable running first over 
			 * 0,...,N_parameters in the ParameterSet, and then over
			 * N_parameters,...,N_parameters+N_observables in the ObservableSet
       * @return The observable. */
	    BCVariable * GetVariable(unsigned int index) const;

      /**
       * @return The number of parameters of the model. */
      unsigned int GetNVariables() const
	       { return fParameters.Size() + fObservables.Size(); }

      /**
       * @param index The index of the parameter in the parameter set.
       * @return The parameter. */
      BCParameter * GetParameter(int index) const
         { return fParameters.Get(index); }

      /**
       * @param name The name of the parameter in the parameter set.
       * @return The parameter. */
      BCParameter * GetParameter(const char * name) const
         { return fParameters.Get(name); }

      /**
       * @return The number of parameters of the model. */
      unsigned int GetNParameters() const
         { return fParameters.Size(); }

      /**
       * @return The number of fixed parameters. */
      unsigned int GetNFixedParameters();

      /**
       * @return The number of free parameters. */
      unsigned int GetNFreeParameters();

      /**
       * @param index The index of the observable in the observable set.
       * @return The user-defined observable. */
      BCObservable * GetObservable(int index) const
         { return fObservables.Get(index); }

      /**
       * @param name The name of the observable in the observable set.
       * @return The user-defined observable. */
      BCObservable * GetObservable(const char * name) const
         { return fObservables.Get(name); }

    	/**
    	 * @return The number of user-defined observables. */
    	unsigned GetNObservables() const
    	   { return fObservables.Size(); }

      /**
       * Returns the set of values of the parameters at the modes of the
       * marginalized posterior pdfs.
       * @return best fit parameters */
      const std::vector<double> & GetBestFitParametersMarginalized() const;

	    /**
			 * Returns the value of a parameter (defined by index) at
			 * the maximum found by the Markov chain
			 * @param index index of the parameter.
			 * @return value of the parameter or -1e+111 on error or center of the range if mode finding not yer run */
	    virtual double GetBestFitParameter(unsigned index) const;

	    /**
			 * @return The log of the value at the mode. */
	    virtual double GetLogMaximum()
	       { return fMCMCLogMaximum; };

	    /**
	     * Returns the value of the user-defined observable at the 
	     * the maximum found by the Markov chain
	     * @param index of the user-defined observable
	     * @return best fit value of user-defined observable. */
	    double GetBestFitObservable(unsigned index) const
	       { return GetBestFitParameter(index+fParameters.Size()); }

      /**
			 * Returns the set of values of the parameters at the maximum
			 * found by the Markov chain
			 * @return best fit parameters */
	    virtual const std::vector<double> & GetBestFitParameters() const
	       { return fMCMCBestFitParameters; }

      /** @} */
      /** \name Setters */
      /** @{ */

      /**
       * Sets the name of the engine.
       * @param name Name of the engine */
	    void SetName(const char * name);

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
       * Sets the number of Markov chains which are run in parallel. */
      void MCMCSetNChains(unsigned n);

      /**
       * Sets the lag of the Markov chains */
      void MCMCSetNLag(unsigned n)
         { fMCMCNLag = n; }

      /**
       * Sets the maximum number of iterations in the pre-run. */
      void MCMCSetNIterationsMax(unsigned n)
         { fMCMCNIterationsMax = n; }

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
      void MCMCSetNIterationsUpdate(unsigned n)
         { fMCMCNIterationsUpdate = n; }

      /**
       * Sets the number of iterations in the pre-run after which the
       * convergence data is cleared.
	   * @param n The number of iterations.*/
      void MCMCSetNIterationsUpdateClear(unsigned n)
         { fMCMCNIterationsUpdateClear = n; }

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
       * Sets flag to write user-defined observables to file during pre run. */
      void MCMCSetWritePreRunObservablesToFile(bool flag)
         { fMCMCFlagWritePreRunObservablesToFile = flag; }

      /**
       * Sets the initial positions for all chains.
       * @param x0s initial positions for all chains. */
      void MCMCSetInitialPositions(const std::vector<double> &x0s);

      /**
       * Sets the initial positions for all chains.
       * @param x0s initial positions for all chains. */
      void MCMCSetInitialPositions(const std::vector< std::vector<double> > & x0s);

      /**
       * Sets flag which defines initial position.  */
      void MCMCSetFlagInitialPosition(int flag)
         { fMCMCFlagInitialPosition = flag; }

      /**
       * Sets the flag which controls the sequence parameters during the
       * running of the MCMC.  */
      void MCMCSetFlagOrderParameters(bool flag)
         { fMCMCFlagOrderParameters = flag; }

      /** Sets the flag for all parameters to either fill histograms or not. */
      void MCMCSetFlagFillHistograms(bool flag);

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

      /** Use strict or relaxed rule for Gelman/Rubin R-value */
      void MCMCSetRValueStrict(bool strict=true)
      { fMCMCRValueUseStrict = strict; }

      /**
       * Sets the tree containing the Markov chains. */
      void MCMCSetMarkovChainTrees(const std::vector<TTree *> & trees);

      /**
       * Initialize trees containing the Markov chains. */
      void MCMCInitializeMarkovChainTrees();

      /**
       * Set the precision for the MCMC run. */
      void MCMCSetPrecision(BCEngineMCMC::Precision precision);

      /**
       * Copy precision for the MCMC run from other model. */
      void MCMCSetPrecision(const BCEngineMCMC * other);
	
#if 0
      /**
       * Set the range of a parameter
       * @param index The parameter index
       * @param parmin The parameter minimum
       * @param parmax The parameter maximum
       * @return An error code. */
      void SetParameterRange(unsigned int index, double parmin, double parmax);
#endif
      /**
       * Set the number of bins for the marginalized distribution of all parameters
       * @param nbins Number of bins */
      void SetNbins(unsigned int nbins);

      /** @} */
      /** \name Error propagation*/
      /** @{ */

      /** @} */

      /**
       * Flag for writing Markov chain to ROOT file (true) or not (false) */
      void WriteMarkovChain(bool flag)
      {
         fMCMCFlagWriteChainToFile = flag;
         fMCMCFlagWritePreRunToFile = flag;
      }

#if 0
      /**
       * Sets the initial position for the Markov chain */
      void SetMarkovChainInitialPosition(const std::vector<double> & position)
      { fXmetro0 = position; }

      /**
       * Sets the step size for Markov chains */
      void SetMarkovChainStepSize(double stepsize)
      { fMarkovChainStepSize = stepsize; }

      /**
       * Sets the number of iterations in the markov chain */
      void SetMarkovChainNIterations(int niterations)
      { fMarkovChainNIterations = niterations;
      fMarkovChainAutoN = false; }

      /**
       * Sets a flag for automatically calculating the number of iterations */
      void SetMarkovChainAutoN(bool flag)
      { fMarkovChainAutoN = flag; }
#endif

      /** @} */
      /** \name Miscellaneous methods */
      /** @{ */

      /**
       * Prints a summary on the screen. */
	    virtual void PrintSummary();

      /**
       * Prints a summary to a file. */
      virtual void PrintResults(const char * file);

	    /**
			 * Print parameters
			 * @P is a vector of the parameter values to be printed
			 * @output is a pointer to the output function to be used,
			 * which defaults to BCLog::OutSummary */
	    void PrintParameters(std::vector<double> const & P, void (*output)(const char *) = BCLog::OutSummary);

	    /**
			 *  Print all 1D marginalizations, each to its own file */
	    int PrintAllMarginalized1D(const char * filebase);

	    /**
			 *  Print all 2D marginalizations, each to its own file */
	    int PrintAllMarginalized2D(const char * filebase);

	    /**
			 *  Print all 1D marginalizations, each to its own file */
	    int PrintAllMarginalized(const char * file, std::string options1d="", std::string options2d="", unsigned int hdiv=1, unsigned int ndiv=1);
	
			/**
			 * Print a summary plot for the parameters and user-defined observables.
			 * @par npar Number of parameters per page, print all on one page if set to zero or negative
			 * @return An error flag. */
			int PrintParameterPlot(const char * filename = "parameters.pdf", int npar=10, double interval_content=68e-2, std::vector<double> quantile_vals=std::vector<double>(0));
			
			/**
			 * Print a summary plot for the parameters in the range provided
			 * @par i0 Index of first parameter to print.
			 * @par npar Number of parameters to print, set to 0 to print all.
			 * @return An error flag. */
			int PrintParameterPlot(unsigned i0, unsigned npar=0, const char * filename = "parameters.pdf", double interval_content=68e-2, std::vector<double> quantile_vals=std::vector<double>(0));
			
			/**
			 * Print a correlation matrix for the parameters.
			 * @return An error flag. */
			int PrintCorrelationMatrix(const char * filename = "matrix.pdf");
			
	    /**
			 * Print a correlation plot for the parameters.
			 * @return An error flag. */
	    int PrintCorrelationPlot(const char * filename = "correlation.pdf");
			
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
      virtual int AddParameter(const char * name, double min, double max, const char * latexname = "");

      /**
       * Adds a parameter to the model.
       * @param parameter A model parameter
       * @see AddParameter(const char * name, double lowerlimit, double upperlimit, const char * latexname); */
      virtual int AddParameter(BCParameter* parameter);

      /**
       * Adds a user-calculated observable.
       * @param min minimum value of the observable
       * @param max maximum value of the observable
       * @param latexname Optional latexname used for plotting
       * @return number of observables after adding */
	    virtual int AddObservable(const char * name, double min, double max, ObservableFunction fn, const char * latexname = "");

      /**
       * Adds a user-calculated observable to the model.
       * @param observable A user-calculated observable
       * @see AddObservable(const char * name, double lowerlimit, double upperlimit, ObservableFunction * fn, const char * latexname); */
      virtual int AddObservable(BCObservable* observable);

      /**
       * Calculates user-defined observables */
  	  virtual void CalculateObservables();

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
       * Updates statistics: find new maximum */
      void MCMCInChainCheckMaximum();

      /**
       * Updates statistics:  */
      void MCMCInChainUpdateStatistics();

      /**
       * Updates statistics: fill marginalized distributions */
      void MCMCInChainFillHistograms();

      /**
       * Updates statistics: check convergence */
      void MCMCInChainTestConvergenceAllChains();

      /**
       * Updates statistics: write chains to file */
      void MCMCInChainWriteChains();

      /**
       * Needs to be overloaded in the derived class.
       * @return natural logarithm of the function to map with MCMC */
      virtual double LogEval(const std::vector<double> & parameters);

      /**
       * Runs Metropolis algorithm. */
      int MCMCMetropolis();

      /**
       * Runs a pre run for the Metropolis algorithm. */
      int MCMCMetropolisPreRun();

      /**
       * Resets the run statistics. */
      void MCMCResetRunStatistics();

      /**
       * Initializes Markov chains. */
      void MCMCInitializeMarkovChains();

      /**
       * Initializes the engine.
       * @return An error code */
      int MCMCInitialize();

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
      virtual void MCMCIterationInterface();

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
	    virtual void PrintSummaryToStream(std::ofstream & ofi);

	    /**
			 * Print best fit to stream. */
	    virtual void PrintBestFitToStream(std::ofstream & ofi);

	    /**
			 * Print marginalization to stream. */
	    virtual void PrintMarginalizationToStream(std::ofstream & ofi);

      /**
       * Name of the engine. */
      std::string fName;

      /**
       * Safe name of the engine for use in naming ROOT objects. */
      std::string fSafeName;

      /**
       * @return The name of the engine with spaces removed. */
      const std::string & GetSafeName() const
         { return fSafeName; }

    	/**
    	 * return appropriate update interval */
    	unsigned  UpdateFrequency(unsigned N);
    
      /**
       * Parameter settings */
      BCParameterSet fParameters;

	    /**
			 * User-calculated Observables Set */
	    BCObservableSet fObservables;

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
      unsigned fMCMCNIterationsUpdate;

      /**
       * Number of iterations for clearing data for updating scale factors */
      unsigned fMCMCNIterationsUpdateClear;

      /**
       * Number of iterations needed for all chains to convergence
       * simultaneously */
      int fMCMCNIterationsConvergenceGlobal;

      /**
       * Maximum number of iterations for a Markov chain prerun */
      unsigned fMCMCNIterationsMax;

      /**
       * Number of iterations for a Markov chain run */
      unsigned fMCMCNIterationsRun;

      /**
       * Minimum number of iterations for the pre-run */
      unsigned fMCMCNIterationsPreRunMin;

      /**
       * Number of accepted trials for each parameter of each chain. */
	    std::vector<std::vector<int> >fMCMCNTrialsTrue;

      /**
       * Number of trials */
      unsigned fMCMCNTrials;

      /**
       * Flag to write Markov chains to file */
      bool fMCMCFlagWriteChainToFile;

      /**
       * Flag to write pre run to file */
      bool fMCMCFlagWritePreRunToFile;

      /**
       * Flag to write user-defined observable to file during pre run. */
      bool fMCMCFlagWritePreRunObservablesToFile;

    	/**
    	 * Lower limit for scale factors */
    	double fMCMCScaleFactorLowerLimit;
    
    	/**
    	 * Upper limit for scale factors */
    	double fMCMCScaleFactorUpperLimit;

      /**
       * Scales the width of the trial functions by a scale factor for
       * each parameter and chain */
	    std::vector<std::vector<double> > fMCMCTrialFunctionScaleFactor;

      /**
       * Start values of the scale factors for the trial functions. */
      std::vector<double> fMCMCTrialFunctionScaleFactorStart;

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
       * The efficiencies for all parameters and chains. */
	    std::vector<std::vector<double> > fMCMCEfficiencies;

      /**
       * The minimum required efficiency for MCMC */
      double fMCMCEfficiencyMin;

      /**
       * The maximum allowed efficiency for MCMC */
      double fMCMCEfficiencyMax;

      /**
       * Variable which defines the initial position. 0 (default) center
       * of the allowed region, (1) random initial position (2)
       * pre-defined intial position. */
      int fMCMCFlagInitialPosition;

      /**
       * Flag which controls the sequence parameters during the running
       * of the MCMC. */
      bool fMCMCFlagOrderParameters;

      /**
       * The phase of the run.
       * 1: pre-run, 2: main run.
       */
      int fMCMCPhase;

      /**
       * The current points of each Markov chain. */
	    std::vector<std::vector<double> > fMCMCx;

      /**
       * The current values of the user-defined observables for each
       * Markov chain. */
	    std::vector<std::vector<double> > fMCMCObservables;

      /**
       * The maximum points of each Markov chain. */
	    std::vector<std::vector<double> > fMCMCxMax;

      /**
       * The mean of all parameters of each Markov chain. */
	    std::vector<std::vector<double> > fMCMCxMean;

      /**
       * The variance of all parameters of each Markov chain. */
	    std::vector<std::vector<double> > fMCMCxVar;

      /**
       * The log of the probability of the current points of each Markov
       * chain. The length of the vectors is fMCMCNChains. */
      std::vector<double> fMCMCprob;

      /**
       * The maximum (log) probability of each Markov chain. The length of
       * the vector is fMCMCNChains. */
      std::vector<double> fMCMCprobMax;

      /**
       * The mean of all log prob values of each Markov chain. The
       * length of the vector is equal to fMCMCNChains. */
      std::vector<double> fMCMCprobMean;

      /**
       * The variance of all log prob values of each Markov chain. The
       * length of the vector is equal to fMCMCNChains. */
      std::vector<double> fMCMCprobVar;

      /**
       * flag: use exactly the R-value definition of Gelman/Rubin (R_strict)
       * or a relaxed, simplified version (R_relax) [default].
       * Note that R_relax <= R_strict, and in some cases even
       * R_relax < 1, but we always have R_strict >= 1.
       */
      bool fMCMCRValueUseStrict;

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
      std::vector<TH1D *>               fH1Marginalized;
	    std::vector<std::vector<TH2D *> > fH2Marginalized;

      /**
       * The trees containing the Markov chains. The length of the vector
       * is fMCMCNChains. */
      std::vector<TTree *> fMCMCTrees;

      /**
       * A vector of best fit parameters found by MCMC */
      std::vector<double> fMCMCBestFitParameters;

      /**
       * The function value at the mode on the @em log scale */
      double fMCMCLogMaximum;

      /**
       * A vector of best fit parameters estimated from the marginalized probability */
      std::vector<double> fMarginalModes;
};

// ---------------------------------------------------------

#endif
