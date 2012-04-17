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

/**
 * Copyright (C) 2008-2012, Daniel Kollar and Kevin Kroeninger.
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

// ---------------------------------------------------------

#include <vector>

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
     enum Precision{ kLow, kMedium, kHigh, kVeryHigh };

      /** @} */
      /** \name Constructors and destructors */
      /** @{ */

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
       * @return number of parameters of the Markov chain */
      int MCMCGetNParameters()
         { return fMCMCNParameters; }

      /**
       * @return number of Markov chains */
      int MCMCGetNChains()
         { return fMCMCNChains; }

      /**
       * @return lag of the Markov chains */
      int MCMCGetNLag()
         { return fMCMCNLag; }

      /**
       * @return number of iterations */
      std::vector<int> MCMCGetNIterations()
         { return fMCMCNIterations; }

      /**
       * @return current iterations */
      int MCMCGetCurrentIteration()
         { return fMCMCCurrentIteration; }

      /**
       * @return current chain index */
      int MCMCGetCurrentChain()
         { return fMCMCCurrentChain; }

      /**
       * @return number of iterations needed for all chains to
       * converge simultaneously */
      int MCMCGetNIterationsConvergenceGlobal()
         { return fMCMCNIterationsConvergenceGlobal; }

      /**
       * @return flag if converged or not */
      bool MCMCGetFlagConvergenceGlobal()
         { return fMCMCFlagConvergenceGlobal; }

      /**
       * @return maximum number of iterations for a Markov chain */
      int MCMCGetNIterationsMax()
         { return fMCMCNIterationsMax; }

      /**
       * @return number of iterations for a Markov chain */
      int MCMCGetNIterationsRun()
         { return fMCMCNIterationsRun; }

      /**
       * @return minimum number of pre-run iterations for a Markov chain */
      int MCMCGetNIterationsPreRunMin()
         { return fMCMCNIterationsPreRunMin; }

      /**
       * @return number of iterations after statistics update. */
      int MCMCGetNIterationsUpdate()
         { return fMCMCNIterationsUpdate; }

      /**
       * @return maximum number of iterations after statistics update. */
      int MCMCGetNIterationsUpdateMax()
         { return fMCMCNIterationsUpdateMax; }

      /**
       * @returns number of accepted trials for each chain */
      std::vector<int> MCMCGetNTrialsTrue()
         { return fMCMCNTrialsTrue; }

      /**
       * @returns number of not-accepted trials for each chain */
      std::vector<int> MCMCGetNTrialsFalse()
         { return fMCMCNTrialsFalse; }

      /**
       * @return mean value of the probability for each chain up to
       * the current iteration  */
      std::vector<double> MCMCGetprobMean()
         { return fMCMCprobMean; }

      /**
       * @return mean value of the probability for each chain up to
       * the current iteration */
      std::vector<double> MCMCGetVariance()
         { return fMCMCprobVar; }

      /**
       * @return scale factor for all parameters and chains */
      std::vector<double> MCMCGetTrialFunctionScaleFactor()
         { return fMCMCTrialFunctionScaleFactor; }

      /**
       * @return scale factor for all parameters and achain.
       * @param ichain chain index */
      std::vector<double> MCMCGetTrialFunctionScaleFactor(int ichain);

      /**
       * @return scale factor for a parameter and a chain.
       * @param ichain chain index
       * @param ipar parameter index */
      double MCMCGetTrialFunctionScaleFactor(int ichain, int ipar);

      /**
       * @return current point of each Markov chain */
      std::vector<double> MCMCGetx()
         { return fMCMCx; }

      /**
       * @param ichain index of the Markov chain
       * @return current point of the Markov chain */
      std::vector<double> MCMCGetx(int ichain);

      /**
       * @param ichain chain index
       * @param ipar parameter index
       * @return parameter of the Markov chain */
      double MCMCGetx(int ichain, int ipar);

      /**
       * @return log of the probability of the current points of each Markov chain */
      std::vector<double> MCMCGetLogProbx()
         { return fMCMCprob; }

      /**
       * @return log of the probability of the current points of the Markov chain.
       * @param ichain chain index */
      double MCMCGetLogProbx(int ichain);

      /**
       * @return pointer to the phase of a run. */
      int MCMCGetPhase()
         { return fMCMCPhase; }

      /**
       * @return pointer to the cycle of a pre-run. */
      int MCMCGetCycle()
         { return fMCMCCycle; }

      /**
       * @return maximum points of each Markov chain */
      std::vector<double> MCMCGetMaximumPoints()
         { return fMCMCxMax; }

      /**
       * @return maximum point of  Markov chain
       * @param i The index of the Markov chain */
      std::vector<double> MCMCGetMaximumPoint(int i);

      /**
       * @return maximum (log) probability of each Markov chain */
      std::vector<double> MCMCGetMaximumLogProb()
         { return fMCMCprobMax; }

      /**
       * @return flag which defined initial position */
      int MCMCGetFlagInitialPosition()
         { return fMCMCFlagInitialPosition; }

      /**
       * @return R-value criterion */
      double MCMCGetRValueCriterion()
         { return fMCMCRValueCriterion; }

      /**
       * @return R-value criterion for parameters */
      double MCMCGetRValueParametersCriterion()
         { return fMCMCRValueParametersCriterion; }

      /**
       * @return R-value */
      double MCMCGetRValue()
         { return fMCMCRValue; }

      /**
       * @return R-value for a parameter
       * @param i parameter index */
      double MCMCGetRValueParameters(int i)
         { return fMCMCRValueParameters.at(i); }

      /** Use strict or relaxed rule for Gelman/Rubin R-value */
      bool MCMCGetRValueStrict()
         { return fMCMCRValueUseStrict; }

      /**
       * @return the flag if MCMC has been performed or not */
      bool MCMCGetFlagRun()
      { return fMCMCFlagRun; }

      /**
       * Rtrieve the tree containing the Markov chain.
       * @param i index of the Markov chain
       * @return pointer to the tree */
      TTree * MCMCGetMarkovChainTree(int i)
         { return fMCMCTrees.at(i); }

      /**
       * Retrieve a histogram of the 1D marginalized distribution of a single parameter.
       * @param i index of the parameter
       * @return pointer to the histogram */
      TH1D * MCMCGetH1Marginalized(int i);

      /**
       * Retrieve a histogram of the 2D marginalized distribution for two parameters.
       * @param i index of the first parameter
       * @param j index of the second parameter
       * @return pointer to the histogram */
      TH2D * MCMCGetH2Marginalized(int i, int j);

      /**
       * Return the random number generator.
       * Any  non-zero seed gives reproducible behavior,e.g.
       * m->MCMCGetTRandom3()->SetSeed(21340) */
      TRandom3 * MCMCGetTRandom3()
         { return fRandom; }

      /** @} */
      /** \name Setters */
      /** @{ */

      /**
       * Set the scale factors for the trial functions
       * @param scale a vector of doubles containing the scale factors */
      void MCMCSetTrialFunctionScaleFactor(std::vector<double> scale)
         { fMCMCTrialFunctionScaleFactorStart = scale; }

      /**
       * Sets the number of Markov chains which are run in parallel. */
      void MCMCSetNChains(int n);

      /**
       * Sets the lag of the Markov chains */
      void MCMCSetNLag(int n)
         { fMCMCNLag = n; }

      /**
       * Sets the maximum number of iterations in the pre-run. */
      void MCMCSetNIterationsMax(int n)
         { fMCMCNIterationsMax = n; }

      /**
       * Sets the number of iterations. */
      void MCMCSetNIterationsRun(int n)
         { fMCMCNIterationsRun = n; }

      /**
       * Sets the minimum number of iterations in the pre-run */
      void MCMCSetNIterationsPreRunMin(int n)
         { fMCMCNIterationsPreRunMin = n; }

      /**
       * Sets the number of iterations in the pre-run after which an
       * update on the statistics (convergence, efficiency, etc.) is done.
       * @param n The number of iterations.*/
      void MCMCSetNIterationsUpdate(int n)
         { fMCMCNIterationsUpdate = n; }

      /**
       * Sets the maximum number of iterations in the pre-run after which an
       * update on the statistics (convergence, efficiency, etc.) is done.
       * If set to 0 no maximum is set.
       * @param n maximum number of iterations. */
      void MCMCSetNIterationsUpdateMax(int n)
         { fMCMCNIterationsUpdateMax = n; }

      /**
       * Sets the minimum efficiency required for a chain. */
      void MCMCSetMinimumEfficiency(double efficiency)
         { fMCMCEfficiencyMin = efficiency; }

      /**
       * Sets the maximum efficiency required for a chain. */
      void MCMCSetMaximumEfficiency(double efficiency)
         { fMCMCEfficiencyMax = efficiency; }

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
      void MCMCSetInitialPositions(std::vector< std::vector<double> > x0s);

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

      /** Sets the flag for a single parameter to either fill histograms or not. */
      void MCMCSetFlagFillHistograms(int index, bool flag);

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
      void MCMCSetMarkovChainTrees(std::vector<TTree *> trees);

      /**
       * Initialize trees containing the Markov chains. */
      void MCMCInitializeMarkovChainTrees();

      /**
       * Sets the histogram with 1D marginalized distributions for parameter.
       * @param i index of the parameter
       * @param h pointer to an existing histogram */
      int SetMarginalized(int index, TH1D * h);

      /**
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

      /**
       * Set the precision for the MCMC run. */ 
      void MCMCSetPrecision(BCEngineMCMC::Precision precision);


      /** @} */
      /** \name Miscellaneous methods */
      /** @{ */

      /**
       * Adds a parameter.
       * @param min minimum value of the parameter
       * @param max maximum value of the parameter
       * @return number of parameters after adding */
      int MCMCAddParameter(double min, double max);

      /**
       * Random walk trial function. The default trial function is a
       * Breit-Wigner. It can be overloaded by the user to set the trial
       * function.
       * @param ichain the chain index
       * @param x point with the dimension fMCMCNParameters */
      virtual void MCMCTrialFunction(int ichain, std::vector<double> &x);

      /**
       * Random walk trial function. The default trial function is a
       * Breit-Wigner. It can be overloaded by the user to set the trial
       * function.
       * @param ichain the chain index
       * @param ipar the parameter index
       * @return the unscaled proposal point */
      virtual double MCMCTrialFunctionSingle(int ichain, int ipar);

      /**
       * Returns a trial point for the Metropolis algorithm.
       * @param chain chain index
       * @param x proposal point
       * @return flag indicating whether the new point lies within the allowed range */
      bool MCMCGetProposalPointMetropolis(int chain, std::vector<double> &x);

      /**
       * Returns a trial point for the Metropolis algorithm.
       * @param chain chain index
       * @param x proposal point
       * @return flag indicating whether the new point lies within the allowed range */
      bool MCMCGetProposalPointMetropolis(int chain, int parameter, std::vector<double> &x);

      /**
       * Generates a new point using the Metropolis algorithm.
       * @param chain chain index */
      bool MCMCGetNewPointMetropolis(int chain = 0);
      bool MCMCGetNewPointMetropolis(int chain, int parameter);

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
       * @return An error code */
      int MCMCResetResults();

      /**
       * Interface allowing to execute arbitrary code for each iteration
       * of the MCMC. The frequency of calling this method is influenced
       * by the setup of the Lag and whether or not the MCMC is run with
       * ordered parameters. This method needs to be overloaded in the derived
       * class. */
      virtual void MCMCIterationInterface()
         {}

      /**
       * Interface allowing to execute arbitrary code for each new point
       * of the MCMC. This method needs to be overloaded in the derived
       * class
       * @param point point that was generated and checked
       * @param ichain index of the chain
       * @param accepted flag whether or not the point was accepted for the chain
       */
      virtual void MCMCCurrentPointInterface(std::vector<double> & point, int ichain, bool accepted)
         {
            // suppress warnings for unused parameters
            // with optimization, no code should be generated
            (void)point;
            (void)ichain;
            (void)accepted;
         }

      /** @} */

   private:

      /**
       * Defines a type of a pointer to a member function. */
      typedef bool (BCEngineMCMC::*MCMCPointerToGetProposalPoint) (int chain, std::vector<double> xnew, std::vector<double> xold) const;

      /**
       * Pointer to a member function */
      MCMCPointerToGetProposalPoint fMCMCPointerToGetProposalPoint;

   protected:

      /**
       * Number of parameters */
      int fMCMCNParameters;

      /**
       * Parameter boundaries */
      std::vector<double> fMCMCBoundaryMin;
      std::vector<double> fMCMCBoundaryMax;

      /**
       * Parameter flags for marginalization */
      std::vector<bool> fMCMCFlagsFillHistograms;

      /**
       * Number of Markov chains ran in parallel */
      int fMCMCNChains;

      /**
       * The lag for the Markov Chain */
      int fMCMCNLag;

      /**
       * Number of total iterations of the Markov chains. The length of
       * the vector is equal to fMCMCNChains. */
      std::vector<int> fMCMCNIterations;

      /**
       * The current iteration number. If not called within the running
       * of the algorithm, return -1. */
      int fMCMCCurrentIteration;

      /**
       * The current chain index. If not called within the running of the
       * algorithm, return -1. */
      int fMCMCCurrentChain;

      /**
       * Number of iterations for updating scale factors */
      int fMCMCNIterationsUpdate;

      /**
       * Maximum number of iterations for updating scale factors */
      int fMCMCNIterationsUpdateMax;

      /**
       * Number of iterations needed for all chains to convergence
       * simulaneously */
      int fMCMCNIterationsConvergenceGlobal;

      /**
       * Flag for convergence */
      bool fMCMCFlagConvergenceGlobal;

      /**
       * Maximum number of iterations for a Markov chain prerun */
      int fMCMCNIterationsMax;

      /**
       * Number of iterations for a Markov chain run */
      int fMCMCNIterationsRun;

      /**
       * Minimum number of iterations for the pre-run */
      int fMCMCNIterationsPreRunMin;

      /**
       * Number of accepted trials for each chain. The length of the
       * vector is equal to fMCMCNChains * fMCMCNParameters.  */
      std::vector<int> fMCMCNTrialsTrue;

      /**
       * Number of not accepted trials for each chain. The length of the
       * vector is equal to fMCMCNChains * fMCMCNParameters.  */
      std::vector<int> fMCMCNTrialsFalse;

      /**
       * Flag to write Markov chains to file */
      bool fMCMCFlagWriteChainToFile;

      /**
       * Flag to write pre run to file */
      bool fMCMCFlagWritePreRunToFile;

      /**
       * Scales the width of the trial functions by a scale factor for
       * each parameter and chain */
      std::vector<double> fMCMCTrialFunctionScaleFactor;


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
      std::vector<double> fMCMCInitialPosition;

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
       * Flag which controls fill histograms during main run. */
      bool fMCMCFlagFillHistograms;

      /**
       * The phase of the run.
       * 1: pre-run, 2: main run.
       */
      int fMCMCPhase;

      /**
       * The cycle of the pre-run
       */
      int fMCMCCycle;

      /**
       * The current points of each Markov chain. The length of the
       * vectors is equal to fMCMCNChains * fMCMCNParameters. First, the
       * values of the first Markov chain are saved, then those of the
       * second and so on. */
      std::vector<double> fMCMCx;

      /**
       * The maximum points of each Markov chain. The length of the vector
       * is fMCMCNChains * fMCMCNParameters. First, the values of the
       * first Markov chain are saved, then those of the second and so on. */
      std::vector<double> fMCMCxMax;

      /**
       * The mean of all parameters of each Markov chain. The length of
       * the vector is equal to fMCMCNChains * fMCMCNParameters. */
      std::vector<double> fMCMCxMean;

      /**
       * The variance of all parameters of each Markov chain. The length
       * of the vector is equal to fMCMCNChains * fMCMCNParameters. */
      std::vector<double> fMCMCxVar;

      /**
       * A temporary vector for a single Markov chain */
      std::vector<double> fMCMCxLocal;

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
       * Number of bins per dimension for the marginalized distributions. */
      std::vector<int> fMCMCH1NBins;

      /**
       * An array of marginalized distributions */
      std::vector<TH1D *> fMCMCH1Marginalized;
      std::vector<TH2D *> fMCMCH2Marginalized;

      /**
       * The trees containing the Markov chains. The length of the vector
       * is fMCMCNChains. */
      std::vector<TTree *> fMCMCTrees;
};

// ---------------------------------------------------------

#endif
