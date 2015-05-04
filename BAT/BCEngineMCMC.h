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

#include "BCParameter.h"
#include "BCParameterSet.h"
#include "BCObservable.h"
#include "BCObservableSet.h"
#include "BCLog.h"

#include <vector>
#include <limits>
#include <utility>

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

class BCEngineMCMC
{

public:

    /** \name Enumerators  */
    /** @{ */

    /** An enumerator for the status of a test.
     * See BCEngineMCMC::SetPrecision for full details*/
    enum Precision {
        kLow,
        kQuick,
        kMedium,
        kHigh,
        kVeryHigh
    };

    /** An enumerator for the phase of the Markov chain.
     * Negative values are in pre-run, 0 is unset, positive is main run. */
    enum MCMCPhase {
        kMCMCPreRun     = -1, ///< In pre-run
        kMCMCUnsetPhase =  0, ///< Unest
        kMCMCMainRun    = +1  ///< In main run
    };

    /** An enumerator for markov-chain position initialization. */
    enum MCMCInitialPosition {
        kMCMCInitContinue         = -1, ///< continue markov chain from last value (not yet available)
        kMCMCInitCenter           =  0, ///< select centers of parameter ranges
        kMCMCInitRandomUniform    =  1, ///< randomly distribute uniformly over parameter ranges
        kMCMCInitUserDefined      =  2, ///< initialize to user-provided points
        kMCMCInitRandomPrior      =  3  ///< randomly distribute according to factorized priors
    };

    /** @} */
    /** \name Structs */
    /** @{ */

    /** A struct for holding statistical information about samples. */
    struct MCMCStatistics {

        /** Constructor.
         * @param n_par number of parameters to calculate statistics for.
         * @param n_obs number of observables to calculate statistics for (sans efficiencies). */
        MCMCStatistics(unsigned n_par = 0, unsigned n_obs = 0);

        /** Copy constructor. */
        MCMCStatistics(const BCEngineMCMC::MCMCStatistics& other);

        unsigned n_samples;					                  ///< number of samples used to calculate statistics
        std::vector<double> mean;		                  ///< means of all variables
        std::vector<double> variance;                 ///< variances of all variables
        std::vector<std::vector<double> > covariance; ///< covariances of all pairs of variables
        std::vector<double> minimum;									///< minimum value of variables
        std::vector<double> maximum;									///< maximum value of variables
        double probability_mean;											///< mean of probability
        double probability_variance;									///< variance of probability
        std::vector<double> mode;											///< mode of variables
        double probability_at_mode;                   ///< mode of probability
        unsigned n_samples_efficiency;								///< number of samples used to calculate efficiencies
        std::vector<double> efficiency;								///< efficiencies for each parameter (NB: not stored for observables)

        /** clear all members.
         * @param clear_mode Flag for clearing information about mode*/
        void Clear(bool clear_mode = true, bool clear_effieciency = true);

        /** init all members
         * @param n_par number of parameters
         * @param n_obs number of observables. */
        void Init(unsigned n_par, unsigned n_obs);

        /** reset all members
         * @param reset_mode flag for resetting information about mode.
         * @param reset_efficiency flag for resetting information about efficiencies. */
        void Reset(bool reset_mode = true, bool reset_efficiency = true);

        /** reset efficiencies */
        void ResetEfficiencies();

        /** assignment operator. */
        MCMCStatistics& operator  = (const MCMCStatistics& rhs);
        /** addition assignment operator. */
        MCMCStatistics& operator += (const MCMCStatistics& rhs);

        /** update statistics given new values.
         * @param par Current parameter values.
         * @param obs Current user-define observables values. */
        void Update(double prob, const std::vector<double>& par, const std::vector<double>& obs);
    };

    /** @} */
    /** \name Constructors and destructors */
    /** @{ */

    /**
     * Default constructor. */
    BCEngineMCMC(std::string name = "model");

    /**
     * Copy constructor. */
    BCEngineMCMC(const BCEngineMCMC& enginemcmc);

    /**
     * Read in MCMC constructor.
     * @param filename Path of file holding model.
     * @param name Name of model (file should contain TTree's [name]_mcmc and [name]_parameters.\n
     * if empty string is given, properly matching TTrees are searched for in the file.
     * @param loadObservables Flag for whether to load observables for file (true; default) or to let user respecify observables.*/
    BCEngineMCMC(std::string filename, std::string name, bool loadObservables = true);

    /**
     * Destructor. */
    virtual ~BCEngineMCMC();

    /** @} */
    /** \name Assignment operators */
    /** @{ */

    /**
     * Defaut assignment operator */
    BCEngineMCMC& operator = (const BCEngineMCMC& engineMCMC)
    { Copy(engineMCMC); return *this; }

    /** @} */
    /** \name Getters */
    /** @{ */

    /**
     * @return The name of the engine. */
    const std::string& GetName() const
    { return fName; }

    /**
     * @return The name of the engine with spaces removed. */
    const std::string& GetSafeName() const
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
    const std::vector<unsigned>& MCMCGetNIterations() const
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
     * @return number of iterations between scale adjustments and convergence checking during pre-run. */
    unsigned MCMCGetNIterationsPreRunCheck() const
    { return fMCMCNIterationsPreRunCheck; }

    /**
     * @return number of iterations between clearings of statistics for convergence checking. */
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
    const std::vector<std::vector<double> >& MCMCGetTrialFunctionScaleFactor() const
    { return fMCMCTrialFunctionScaleFactor; }

    /**
     * @return scale factor for all parameters of a chain.
     * @param c chain index */
    const std::vector<double>& MCMCGetTrialFunctionScaleFactor(unsigned c) const
    { return fMCMCTrialFunctionScaleFactor.at(c); }

    /**
     * @return scale factor for a parameter and a chain.
     * @param c chain index
     * @param p parameter index */
    double MCMCGetTrialFunctionScaleFactor(unsigned c, unsigned p) const
    { return fMCMCTrialFunctionScaleFactor.at(c).at(p); }

    /**
     * @return current point of each Markov chain */
    const std::vector<std::vector<double> >& MCMCGetx() const
    { return fMCMCx; }

    /**
     * @param c index of the Markov chain
     * @return current point of the Markov chain */
    const std::vector<double>& MCMCGetx(unsigned c) const
    { return fMCMCx.at(c); }

    /**
     * @param c chain index
     * @param p parameter index
     * @return parameter of the Markov chain */
    double MCMCGetx(unsigned c, unsigned p) const
    { return fMCMCx.at(c).at(p); }

    /**
     * @return log of the probability of the current points of each Markov chain */
    const std::vector<double>& MCMCGetLogProbx() const
    { return fMCMCprob; }

    /**
     * @return log of the probability of the current points of the Markov chain.
     * @param c chain index */
    double MCMCGetLogProbx(unsigned c) const
    { return fMCMCprob.at(c); }

    /**
     * @return pointer to the phase of a run. */
    BCEngineMCMC::MCMCPhase MCMCGetPhase() const
    { return fMCMCPhase; }

    /**
     * @return flag which defined initial position */
    BCEngineMCMC::MCMCInitialPosition MCMCGetFlagInitialPosition() const
    { return fMCMCFlagInitialPosition; }

    /**
     * @return whether to use a multivariate proposal function. */
    bool MCMCGetMultivariateProposalFunction() const
    { return fMCMCMultivariateProposalFunction; }

    /**
     * @return minimum number of updates to multivariate-proposal-function covariance to perform. */
    unsigned MCMCGetMultivariateProposalFunctionUpdatesMinimum() const
    { return fMultivariateProposalFunctionUpdatesMinimum; }

    /**
     * @return multivariate-proposal-function Cholesky decomposition nudge size. */
    double MCMCGetMultivariateProposalFunctionEpsilon() const
    { return fMultivariateProposalFunctionEpsilon; }

    /**
     * @return multivariate-proposal-function scale multiplier. */
    double MCMCGetMultivariateProposalFunctionScaleMultiplier() const
    { return fMultivariateProposalFunctionScaleMultiplier; }

    /**
     * @return R-value criterion for parameters */
    double MCMCGetRValueParametersCriterion() const
    { return fMCMCRValueParametersCriterion; }

    /**
     * @return vector of R values for parameters */
    const std::vector<double>& MCMCGetRValueParameters() const
    { return fMCMCRValueParameters; }

    /**
     * @return R-value for a parameter
     * @param i parameter index */
    double MCMCGetRValueParameters(unsigned p) const
    { return fMCMCRValueParameters.at(p); }

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
    TTree* MCMCGetMarkovChainTree()
    { return fMCMCTree;}

    /**
     * Retrieve the tree containing the parameter information. */
    TTree* MCMCGetParameterTree()
    { return fParameterTree;}

    /**
     * Retrieve output file for MCMC. */
    TFile* MCMCGetOutputFile()
    { return fMCMCOutputFile; }

    /**
     * Get combined statistics for all chains. */
    const BCEngineMCMC::MCMCStatistics& MCMCGetStatistics() const
    { return fMCMCStatistics_AllChains; }

    /**
     * Get MCMC statistics for one chain.
     * @param c Chain to get statistics of. */
    const BCEngineMCMC::MCMCStatistics& MCMCGetStatistics(unsigned c) const
    { return fMCMCStatistics.at(c); }

    /**
     * @return Flag for whether to rescale histogram ranges to fit MCMC reach after pre-run. */
    bool GetRescaleHistogramRangesAfterPreRun() const
    { return fRescaleHistogramRangesAfterPreRun; }

    /**
     * @return Factor by which to enlarge histogram ranges when rescaling to add padding beyond range. */
    double GetHistogramRescalePadding() const
    { return fHistogramRescalePadding; }

    /**
     * @return vector of index pairs for order of printing 2D histograms. */
    virtual std::vector<std::pair<unsigned, unsigned> > GetH2DPrintOrder() const;

    /**
     * @param index Index of histogram of which to check existence
     * @return Whether the marginalized histogram exists. */
    bool MarginalizedHistogramExists(unsigned index) const
    { return index < fH1Marginalized.size() and fH1Marginalized[index]; }

    /**
     * @param index1 X Index of histogram of which to check existence.
     * @param index2 Y Index of histogram of which to check existence.
     * @return Whether the marginalized histogram exists. */
    bool MarginalizedHistogramExists(unsigned index1, unsigned index2) const
    { return index1 < fH2Marginalized.size() and index2 < fH2Marginalized[index1].size() and fH2Marginalized[index1][index2]; }

    /**
     * Obtain the individual marginalized distributions
     * with respect to one parameter as a ROOT TH1
     * @note The most efficient method is to access by index.
     * @param name The parameter's name
     * @return 1D marginalized probability */
    TH1* GetMarginalizedHistogram(std::string name) const
    { return GetMarginalizedHistogram(fParameters.Index(name)); }

    /**
     * Obtain the individual marginalized distributions
     * with respect to one parameter as a ROOT TH1
     * @param index The parameter index
     * @return 1D marginalized probability */
    TH1* GetMarginalizedHistogram(unsigned index) const;

    /**
     * Obtain the individual marginalized distributions
     * with respect to two parameters as a ROOT TH2.
     * @note The most efficient method is to access by indices.
     * @param name1 Name of first parameter
     * @param name2 Name of second parameter
     * @return 2D marginalized probability */
    TH2* GetMarginalizedHistogram(std::string name1, std::string name2) const
    { return GetMarginalizedHistogram(fParameters.Index(name1), fParameters.Index(name2)); }

    /**
     * Obtain the individual marginalized distributions
     * with respect to two parameters as a ROOT TH2.
     * @param index1 Index of first parameter
     * @param index2 Index of second parameter
     * @return 2D marginalized probability */
    TH2* GetMarginalizedHistogram(unsigned index1, unsigned index2) const;

    /**
     * Obtain the individual marginalized distributions
     * with respect to one parameter.
     * @note The most efficient method is to access by index.
     * @note Ownership of the returned heap object is conferred to the caller.
     * @param name The parameter's name
     * @return 1D marginalized probability */
    BCH1D* GetMarginalized(std::string name) const
    { return GetMarginalized(fParameters.Index(name)); }

    /**
     * Obtain the individual marginalized distributions
     * with respect to one parameter.
     * @note Ownership of the returned heap object is conferred to the caller.
     * @param index The parameter index
     * @return 1D marginalized probability */
    BCH1D* GetMarginalized(unsigned index) const;

    /**
     * Obtain the individual marginalized distributions
     * with respect to two parameters.
     * @note The most efficient method is to access by indices.
     * @note Ownership of the returned heap object is conferred to the caller.
     * @param name1 Name of first parameter
     * @param name2 Name of second parameter
     * @return 2D marginalized probability */
    BCH2D* GetMarginalized(std::string name1, std::string name2) const
    { return GetMarginalized(fParameters.Index(name1), fParameters.Index(name2)); }

    /**
     * Obtain the individual marginalized distributions
     * with respect to two parameters.
     * @note Ownership of the returned heap object is conferred to the caller.
     * @param index1 Index of first parameter
     * @param index2 Index of second parameter
     * @return 2D marginalized probability */
    BCH2D* GetMarginalized(unsigned index1, unsigned index2) const;

    /**
     * @param observables Whether to check max length of user-defined observable names.
     * @return length of longest parameter name. */
    unsigned GetMaximumParameterNameLength(bool observables = true) const
    { return (observables) ? std::max(fParameters.MaxNameLength(), fObservables.MaxNameLength()) : fParameters.MaxNameLength(); }

    /**
     * @param index The index of the variable running first over
     * 0,...,N_parameters in the ParameterSet, and then over
     * N_parameters,...,N_parameters+N_observables in the ObservableSet
     * @return The variable. */
    BCVariable& GetVariable(unsigned index)
    { return (index < GetNParameters()) ? static_cast<BCVariable&>(fParameters[index]) : static_cast<BCVariable&>(fObservables.At(index - GetNParameters())); }

    /**
     * @param index The index of the observable running first over
     * 0,...,N_parameters in the ParameterSet, and then over
     * N_parameters,...,N_parameters+N_observables in the ObservableSet
     * @return The observable. */
    const BCVariable& GetVariable(unsigned index) const
    { return (index < GetNParameters()) ? static_cast<const BCVariable&>(fParameters[index]) : static_cast<const BCVariable&>(fObservables.At(index - GetNParameters())); }

    /**
     * @return The number of parameters of the model. */
    unsigned GetNVariables() const
    { return fParameters.Size() + fObservables.Size(); }

    /**
     * @return Parameter set. */
    BCParameterSet& GetParameters()
    { return fParameters; }

    /**
     * @return Parameter set. */
    const BCParameterSet& GetParameters() const
    { return fParameters; }

    /**
     * @deprecated Instead call GetParameters().At(index)
     * @param index The index of the parameter in the parameter set.
     * @return The parameter. */
    BCParameter& GetParameter(unsigned index)
    { return fParameters.At(index); }

    /**
     * @deprecated Instead call GetParameters().At(index)
     * @param index The index of the parameter in the parameter set.
     * @return The parameter. */
    const BCParameter& GetParameter(unsigned index) const
    { return fParameters.At(index); }

    /**
     * @deprecatd Instead call GetParameters().Get(name)
     * @param name The name of the parameter in the parameter set.
     * @return The parameter. */
    BCParameter& GetParameter(const std::string& name)
    { return fParameters.Get(name); }

    /**
     * @deprecated Instead call GetParameters().Get(name)
     * @param name The name of the parameter in the parameter set.
     * @return The parameter. */
    const BCParameter& GetParameter(const std::string& name) const
    { return fParameters.Get(name); }

    /**
     * @return The number of parameters of the model. */
    unsigned GetNParameters() const
    { return fParameters.Size(); }

    /**
     * @deprecated Instead call GetParameters().GetNFixedParameters()
     * @return The number of fixed parameters. */
    unsigned GetNFixedParameters() const
    { return fParameters.GetNFixedParameters(); }

    /**
     * @deprecated Instead call GetParameters().GetNFreeParameters()
     * @return The number of free parameters. */
    unsigned GetNFreeParameters() const
    { return fParameters.GetNFreeParameters(); }

    /**
     * @return Observable set. */
    BCObservableSet& GetObservables()
    { return fObservables; }

    /**
     * @return Observable set. */
    const BCObservableSet& GetObservables() const
    { return fObservables; }

    /**
     * @deprecated Instead call GetObservables().At(index)
     * @param index The index of the observable in the observable set.
     * @return The user-defined observable. */
    BCObservable& GetObservable(unsigned index)
    { return fObservables.At(index); }

    /**
     * @deprecated Instead call GetObservbles().At(index)
     * @param index The index of the observable in the observable set.
     * @return The user-defined observable. */
    const BCObservable& GetObservable(unsigned index) const
    { return fObservables.At(index); }

    /**
     * @deprecated Instead call GetObservables().Get(name)
     * @param name The name of the observable in the observable set.
     * @return The user-defined observable. */
    BCObservable& GetObservable(const std::string& name)
    { return fObservables.Get(name); }

    /**
     * @deprecated Instead call GetObservables().Get(name)
     * @param name The name of the observable in the observable set.
     * @return The user-defined observable. */
    const BCObservable& GetObservable(const std::string& name) const
    { return fObservables.Get(name); }

    /**
     * @return The number of user-defined observables. */
    unsigned GetNObservables() const
    { return fObservables.Size(); }

    /**
     * @return vector of parameter and observable values at global mode. */
    virtual const std::vector<double>& GetGlobalMode() const
    { return fMCMCStatistics_AllChains.mode; }

    /**
     * @return vector of the local modes of parameters and observables
     * @param force_recalculation flag for forcing recalculation of local modes from histograms. */
    const std::vector<double>& GetLocalModes(bool force_recalculation = false);

    /**
     * @return The log of the value at the mode. */
    virtual double GetLogMaximum() const
    { return fMCMCStatistics_AllChains.probability_at_mode; };

    /**
     * @return Flag whether to reuse user-defined observables from MCMC tree when looping through it. */
    bool MCMCGetReuseObservables() const
    { return fMCMCTreeReuseObservables; }

    /**
     * @return BCH1D object that stores drawing options for all BCH1D's. */
    BCH1D* GetBCH1DdrawingOptions()
    { return fBCH1DdrawingOptions; }

    /**
     * @return BCH2D object that stores drawing options for all BCH2D's. */
    BCH2D* GetBCH2DdrawingOptions()
    { return fBCH2DdrawingOptions; }

    /** @} */
    /** \name Setters */
    /** @{ */

    /**
     * Sets the name of the engine.
     * @param name Name of the engine */
    void SetName(std::string name);

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
    void MCMCAutoSetTrialFunctionScaleFactors(bool flag = true)
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
     * Sets the number of iterations between scale adjustments and
     * convergence checks in the pre-run. */
    void MCMCSetNIterationsPreRunCheck(unsigned n)
    { fMCMCNIterationsPreRunCheck = n; }

    /**
     * Sets the number of iterations between clearings of convergence
     * checking stats. */
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
    void MCMCSetInitialPositions(const std::vector<double>& x0s);

    /**
     * Sets the initial positions for all chains.
     * @param x0s initial positions for all chains. */
    void MCMCSetInitialPositions(const std::vector< std::vector<double> >& x0s)
    { fMCMCInitialPosition = x0s; MCMCSetFlagInitialPosition(BCEngineMCMC::kMCMCInitUserDefined); }

    /**
     * Sets flag which defines initial position.  */
    void MCMCSetFlagInitialPosition(BCEngineMCMC::MCMCInitialPosition flag)
    { fMCMCFlagInitialPosition = flag; }

    /**
     * Sets the flag which controls the sequence parameters during the
     * running of the MCMC.  */
    void MCMCSetMultivariateProposalFunction(bool flag)
    { fMCMCMultivariateProposalFunction = flag; }

    /**
     * Set minimum number of updates to multivariate-proposal-function covariance to perform. */
    void MCMCSetMultivariateProposalFunctionUpdatesMinimum(unsigned n)
    { fMultivariateProposalFunctionUpdatesMinimum = n; }

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
    { fParameters.FillHistograms(flag_1d, flag_2d); fObservables.FillHistograms(flag_1d, flag_2d); }

    /** Sets whether to fill particular H2 histogram: par(y) vs. par(x)
     * @param x Index of parameter for horizontal axis.
     * @param y Index of parameter for vertical axis.
     * @param flag Whether to fill histogram. */
    void SetFillHistogramParPar(unsigned x, unsigned y, bool flag = true)
    { SetFillHistogram(x, y, flag); }

    /** Sets whether to fill particular H2 histogram: par(y) vs. par(x)
     * @param x Name of parameter for horizontal axis.
     * @param y Name of parameter for vertical axis.
     * @param flag Whether to fill histogram. */
    void SetFillHistogramParPar(std::string x, std::string y, bool flag = true)
    { SetFillHistogramParPar(fParameters.Index(x), fParameters.Index(y), flag); }

    /** Sets whether to fill particular H2 histogram: obs(y) vs. par(x)
     * @param x Index of parameter for horizontal axis.
     * @param y Index of observable for vertical axis.
     * @param flag Whether to fill histogram. */
    void SetFillHistogramParObs(unsigned x, unsigned y, bool flag = true)
    { SetFillHistogram(x, -(y + 1), flag); }

    /** Sets whether to fill particular H2 histogram: obs(y) vs. par(x)
     * @param x Name of parameter for horizontal axis.
     * @param y Name of observable for vertical axis.
     * @param flag Whether to fill histogram. */
    void SetFillHistogramParObs(std::string x, std::string y, bool flag = true)
    { SetFillHistogramParObs(fParameters.Index(x), fObservables.Index(y), flag); }

    /** Sets whether to fill particular H2 histogram: obs(y) vs. obs(x)
     * @param x Index of observable for horizontal axis.
     * @param y Index of observable for vertical axis.
     * @param flag Whether to fill histogram. */
    void SetFillHistogramObsObs(unsigned x, unsigned y, bool flag = true)
    { SetFillHistogram(-(x + 1), -(y + 1), flag); }

    /** Sets whether to fill particular H2 histogram: obs(y) vs. obs(x)
     * @param x Name of observable for horizontal axis.
     * @param y Name of observable for vertical axis.
     * @param flag Whether to fill histogram. */
    void SetFillHistogramObsObs(std::string x, std::string y, bool flag = true)
    { SetFillHistogramObsObs(fObservables.Index(x), fObservables.Index(y), flag); }

    /** Sets whether to fill particular H2 histogram: par(y) vs. obs(x)
     * @param x Index of observable for horizontal axis.
     * @param y Index of parameter for vertical axis.
     * @param flag Whether to fill histogram. */
    void SetFillHistogramObsPar(unsigned x, unsigned y, bool flag = true)
    { SetFillHistogram(-(x + 1), y, flag); }

    /** Sets whether to fill particular H2 histogram: par(y) vs. obs(x)
     * @param x Name of observable for horizontal axis.
     * @param y Name of parameter for vertical axis.
     * @param flag Whether to fill histogram. */
    void SetFillHistogramObsPar(std::string x, std::string y, bool flag = true)
    { SetFillHistogramObsPar(fObservables.Index(x), fParameters.Index(y), flag); }

    /** Sets the flag if a prerun should be performed or not. */
    void MCMCSetFlagPreRun(bool flag)
    { fMCMCFlagPreRun = flag; }

    /**
     * Sets the parameter R-value criterion for convergence of all chains */
    void MCMCSetRValueParametersCriterion(double r)
    { fMCMCRValueParametersCriterion = r; }

    /** Set flag to correct convergence checking for initial sampling variability. */
    void MCMCSetCorrectRValueForSamplingVariability(bool flag = true)
    { fCorrectRValueForSamplingVariability = flag; }

    /**
     * Set the precision for the MCMC run. */
    void MCMCSetPrecision(BCEngineMCMC::Precision precision);

    /**
     * Copy precision for the MCMC run from other model. */
    void MCMCSetPrecision(const BCEngineMCMC* other);

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

    /** @} */

    /** \name Prior setting functions (all deprecated).
     * \brief The priors are not used for the likelihood calculation by
     * BCEngineMCMC, but are used for initializing the positions of the
     * chains. */
    /** @{ */

    /**
     * @deprecated Instead call: GetParameter(index)->SetPriorConstant()
     * Set constant prior for this parameter
     * @param index the index of the parameter */
    void SetPriorConstant(unsigned index)
    { fParameters.At(index).SetPriorConstant(); }

    /**
     * @deprecated Instead call: GetParameter(name)->SetPriorConstant()
     * Set constant prior for this parameter
     * @param name the name of the parameter */
    void SetPriorConstant(std::string name)
    { SetPriorConstant(fParameters.Index(name)); }

    /**
     * @deprecated Instead call: GetPrarameter(index)->SetPrior(new BCTF1Prior(f))
     * Set prior for a parameter.
     * @param index The parameter index
     * @param f A pointer to a function describing the prior */
    void SetPrior(unsigned index, TF1* const f);

    /**
     * @deprecated Instead call: GetParameter(name)->SetPrior(new BCTF1Prior(f))
     * Set prior for a parameter.
     * @param name The parameter name
     * @param f A pointer to a function describing the prior */
    void SetPrior(std::string name, TF1* const f)
    { return SetPrior(fParameters.Index(name), f); }

    /**
     * @deprecated Instead call: GetParameter(index)->Fix(value)
     * Fixes parameter to value.
     * @param index The parameter index
     * @param value The position of the delta function. */
    void SetPriorDelta(unsigned index, double value)
    { fParameters.At(index).Fix(value); }

    /**
     * @deprecated Instead call: GetParameter(name)->Fix(value)
     * Fixes parameter to value.
     * @param name The parameter name
     * @param value The position of the delta function. */
    void SetPriorDelta(std::string name, double value)
    { SetPriorDelta(fParameters.Index(name), value); }

    /**
     * @deprecated Instead call: GetParameter(index)->SetPrior(new BCGaussianPrior(mean,sigma))
     * Set Gaussian prior for a parameter.
     * @param index The parameter index
     * @param mean The mean of the Gaussian
     * @param sigma The sigma of the Gaussian */
    void SetPriorGauss(unsigned index, double mean, double sigma);

    /**
     * @deprecated Instead call: GetParameter(name)->SetPrior(new BCGaussianPrior(mean,sigma))
     * Set Gaussian prior for a parameter.
     * @param name The parameter name
     * @param mean The mean of the Gaussian
     * @param sigma The sigma of the Gaussian */
    void SetPriorGauss(std::string name, double mean, double sigma)
    { SetPriorGauss(fParameters.Index(name), mean, sigma); }

    /**
     * @deprecated Instead call: GetParameter(index)->SetPrior(new BCSplitGaussianPrior(mode,sigma_below,sigma_above))
     * Set Gaussian prior for a parameter with two different widths.
     * @param index The parameter index
     * @param mode The mode of the Gaussian
     * @param sigma_below Standard deviation below mode.
     * @param sigma_above Standard deviation above mode. */
    void SetPriorGauss(unsigned index, double mode, double sigma_below, double sigma_above);

    /**
     * @deprecated Instead call: GetParameter(name)->SetPrior(new BCSplitGaussianPrior(mode,sigma_below,sigma_above))
     * Set Gaussian prior for a parameter with two different widths.
     * @param name The parameter name
     * @param mode The mode of the Gaussian
     * @param sigma_below Standard deviation below mode.
     * @param sigma_above Standard deviation above mode. */
    void SetPriorGauss(std::string name, double mode, double sigmadown, double sigmaup)
    {	SetPriorGauss(fParameters.Index(name), mode, sigmadown, sigmaup); }

    /**
     * @deprecated Instead call: GetParameter(index)->SetPrior(new BCTH1Prior(h,interpolate))
     * Set prior for a parameter.
     * @param index parameter index
     * @param h pointer to a histogram describing the prior
     * @param interpolate whether or not to use linear interpolation */
    void SetPrior(unsigned index, TH1* const h, bool interpolate = false);

    /**
     * @deprecated Instead call: GetParameter(name)->SetPrior(new BCTH1Prior(h,interpolate))
     * Set prior for a parameter.
     * @param name parameter name
     * @param h pointer to a histogram describing the prior
     * @param interpolate whether or not to use linear interpolation
     * @return success of action. */
    void SetPrior(std::string name, TH1* const h, bool interpolate = false)
    { SetPrior(fParameters.Index(name), h, interpolate); }

    /**
     * @deprecated Instead call: GetParameters().SetPriorConstantAll() */
    void SetPriorConstantAll()
    { fParameters.SetPriorConstantAll(); }

    /** @} */

    /** \name Output functions */
    /** @{ */

    /**
     * Write marginalization histograms to file.
     * @param filename Path to write file to.
     * @param option Options passed to ROOT's TFile::Open. */
    void WriteMarginalizedDistributions(std::string filename, std::string option);

    /**
     * Prints a summary on the screen. */
    virtual void PrintSummary();

    /**
     * Prints a summary to a file.
     * @param file Path to file to print results to.
     * @return Success of action.*/
    virtual bool PrintResults(std::string file) const;

    /**
     * Print parameters
     * @param P vector of the parameter values to be printed
     * @param output pointer to the output function to be used, which defaults to BCLog::OutSummary */
    void PrintParameters(const std::vector<double>& P, void (*output)(std::string) = BCLog::OutSummary) const;

    /**
     * Print all marginalizations.
     * @param filename Path to file to print to
     * @param hdiv Number of columns of plots per page
     * @param vdiv Number of rows of plots per page
     * @return Number of plots printed */
    unsigned PrintAllMarginalized(std::string filename, unsigned hdiv = 1, unsigned vdiv = 1) const;

    /**
     * Print a summary plot for the parameters and user-defined observables.
     * @param filename Path to filename to print to.
     * @param npar Number of parameters per page, print all on one page if set to zero or negative
     * @param interval_content Probability mass to display in smallest X interval band
     * @param quantile_values Vector of quantile values to draw
     * @param rescale_ranges Flag for rescaling to range surveyed by MCMC chains
     * @return Number of pages printed. */
    unsigned PrintParameterPlot(std::string filename, int npar = 10, double interval_content = 68e-2, std::vector<double> quantile_vals = std::vector<double>(0), bool rescale_ranges = true) const;

    /**
     * Draw a summary plot for the parameters in the range provided to current pad
     * @par i0 Index of first parameter to print.
     * @par npar Number of parameters to print, set to 0 to print all.
     * @param interval_content Probability mass to display in smallest X interval band
     * @param quantile_values Vector of quantile values to draw
     * @param rescale_ranges Flag for rescaling to range surveyed by MCMC chains
     * @return Success of action. */
    bool DrawParameterPlot(unsigned i0, unsigned npar = 0, double interval_content = 68e-2, std::vector<double> quantile_vals = std::vector<double>(0), bool rescale_ranges = true) const;

    /**
     * Print a correlation matrix for the parameters.
     * @param filename Path to file to print correlation matrix to.
     * @return Success of action. */
    bool PrintCorrelationMatrix(std::string filename = "matrix.pdf") const;

    /**
     * Print a correlation plot for the parameters.
     * @param filename Path to file to print correlation plot to.
     * @param include_observables Flag for including observables (default: true)
     * @return Success of action. */
    bool PrintCorrelationPlot(std::string filename = "correlation.pdf", bool include_observables = true) const;

    /**
     * Print a LaTeX table of the parameters.
     * @param filename Path to file tp print LaTeX table of parameters to.
     * @return Success of action. */
    bool PrintParameterLatex(std::string filename) const;

    /** @} **/
    /** \name Miscellaneous methods */
    /** @{ */

    /**
     * Create histograms from parameter and observable sets.
     * @param rescale_ranges Rescale axis ranges to range reached by MCMC tree if true */
    virtual void CreateHistograms(bool rescale_ranges = false);

    /**
     * Initialize the trees containing the Markov chains and parameter info.
     * @param replacetree Flag to delete and recreate tree object if already existing.
     * @param replacefile Flag to delete and recreate file object if already existing.*/
    virtual void InitializeMarkovChainTree(bool replacetree = false, bool replacefile = false);

    /**
     * Copy object
     * @param enginemcmc Object to copy from */
    void Copy(const BCEngineMCMC& enginemcmc);

    /**
     * @deprecated Instead use GetParameters().Add(...)
     * Adds a parameter.
     * @param name Name of parameter
     * @param min minimum value of the parameter
     * @param max maximum value of the parameter
     * @param latexname Optional latexname used for plotting
     * @param unitstring Unit string to be printed for parameter.
     * @return Success of action. */
    virtual bool AddParameter(std::string name, double min, double max, std::string latexname = "", std::string unitstring = "")
    { return fParameters.Add(name, min, max, latexname, unitstring); }

    /**
     * @deprecated Instead use GetParameters().Add(parameter)
     * Adds a parameter to the model.
     * @param parameter A model parameter
     * @return Success of action. */
    virtual bool AddParameter(BCParameter& parameter)
    { return fParameters.Add(parameter); }

    /**
     * @deprecated Instead use GetObservables().Add(...)
     * Adds a user-calculated observable.
     * @param name name of observable
     * @param min minimum value of the observable
     * @param max maximum value of the observable
     * @param latexname Optional latexname used for plotting
     * @param unitstring Unit string to be printed for observable.
     * @return Success of action. */
    virtual bool AddObservable(std::string name, double min, double max, std::string latexname = "", std::string unitstring = "")
    { return fObservables.Add(name, min, max, latexname, unitstring); }

    /**
     * @deprecated Instead use GetObservables().Add(obs)
     * Adds a user-calculated observable to the model.
     * @param observable A user-calculated observable
     * @return Success of action. */
    virtual bool AddObservable(BCObservable& obs)
    { return fObservables.Add(obs); }

    /**
     * Evaluates user-defined observables at current state of all
     * chains and stores results in fMCMCObservables*/
    virtual void EvaluateObservables();

    /**
     * Evaluates user-defined observables at current state of chain
     * and stores results in fMCMCObservables
     * @param chain Chain to evaluate observables for. */
    virtual void EvaluateObservables(unsigned chain);

    /**
     * Evaluates user-defined observables. To be overloaded by user
     * to calculate user-observables.
     * @param pars Parameter set to evaluate observables for. */
    virtual void CalculateObservables(const std::vector<double>& /*pars*/)
    {}

    /**
     * Random walk trial function. The default trial function is a Breit-Wigner.
     * It can be overloaded by the user to set the trial function.
     * @param ichain the chain index
     * @param x point with the dimension fMCMCNParameters */
    virtual void MCMCTrialFunction(unsigned ichain, std::vector<double>& x);

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
    bool MCMCGetProposalPointMetropolis(unsigned chain, std::vector<double>& x);

    /**
     * Returns a trial point for the Metropolis algorithm.
     * @param chain chain index
     * @param x proposal point
     * @return flag indicating whether the new point lies within the allowed range */
    bool MCMCGetProposalPointMetropolis(unsigned chain, unsigned parameter, std::vector<double>& x);

    /**
     * Generates a new point using the Metropolis algorithm for all chains.
     * @return Whether proposed point was accepted (true) or previous point was kept (false). */
    bool MCMCGetNewPointMetropolis();

    /** Generates a new point using the Metropolis algorithm for one chain.
     * @param chain chain index
     * @return Whether proposed point was accepted (true) or previous point was kept (false). */
    bool MCMCGetNewPointMetropolis(unsigned chain);

    /** Generates a new point using the Metropolis algorithm for one chain, varying only one parameter's value
     * @param chain chain index
     * @param parameter index of parameter to vary
     * @return Whether proposed point was accepted (true) or previous point was kept (false). */
    bool MCMCGetNewPointMetropolis(unsigned chain, unsigned parameter);

    /**
     * Updates statistics: fill marginalized distributions */
    void MCMCInChainFillHistograms();

    /**
     * Updates statistics: write chains to file */
    void MCMCInChainWriteChains();

    /**
     * Needs to be overloaded in the derived class.
     * @param parameters Parameter set to evaluate at.
     * @return natural logarithm of the function to map with MCMC */
    virtual double LogEval(const std::vector<double>& parameters) = 0;

    /**
     * Runs Metropolis algorithm.
     * @return Success of action. */
    bool MCMCMetropolis();

    /**
     * Runs a pre run for the Metropolis algorithm.
     * @return Success of action. */
    bool MCMCMetropolisPreRun();

    /**
     * Calculate R value of set of batches of samples---represented by
     * their means and variances, all batches containing the same
     * number of samples---according to Brooks & Gelman, "General
     * Methods for Monitoring Convergence of Iterative Simulations," (1988)
     * @param means Vector of means of sample batches.
     * @param variances Vector of variances of samples batches.
     * @param n Number of samples in each batch.
     * @param correctForSamplingVariability Flag to control correcting R value for initial sampling variability.
     * @return R value for set of batches of samples. */
    static double RValue(std::vector<double> means, std::vector<double> variances, unsigned n, bool correctForSamplingVariability = true);

    /**
     * Resets all containers used in MCMC and initializes starting points.
     * @return Success of action. */
    bool MCMCInitialize();

    /**
     * Reset the MCMC variables. */
    virtual void ResetResults();

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
     * @param accepted flag whether or not the point was accepted for the chain */
    virtual void MCMCCurrentPointInterface(std::vector<double>& /*point*/, int /*ichain*/, bool /*accepted*/)
    {}

    /**
     * Load parameters and observables from tree.
     * @param partree Tree holding parameter information.
     * @param loadObservables Flag for whether to also load observables.
     * @return Success of action. */
    virtual bool LoadParametersFromTree(TTree* partree, bool loadObservables = true);

    /**
     * Check parameter tree against model.
     * @param partree Tree of parameters to check against model.
     * @param checkObservables Flag for whether to check observables.
     * @param Whether tree's parameters match those of model. */
    virtual bool ParameterTreeMatchesModel(TTree* partree, bool checkObservables = true);

    /**
     * Load previous MCMC run.
     * @param filename Pathname of file containing model trees.
     * @param mcmcTreeName Name of tree inside file containing MCMC, empty string (default) loads [modelname]_mcmc.
     * @param parameterTreeName Name of tree inside file containing parameter list, empty string (default) loads [modelname]_parameters.
     * @param loadObservables Flag for whether to load observables from parameter list and MCMC trees. */
    virtual bool LoadMCMC(std::string filename, std::string mcmcTreeName = "", std::string parameterTreeName = "", bool loadObservables = true);

    /**
     * Load previous MCMC run.
     * @param mcmcTree Tree containing MCMC samples.
     * @param parTree Tree containing definition of parameters.
     * @param loadObservables Flag for whether to load observables.
     * @return Success of action. */
    virtual bool LoadMCMC(TTree* mcmcTree, TTree* parTree, bool loadObservables = true);

    /**
     * Check tree structure for MCMC tree.
     * @param tree MCMC Tree to check validity of.
     * @param checkObservables Flag for whether to check observables.
     * @return Validity of tree. */
    virtual bool ValidMCMCTree(TTree* tree, bool checkObservables = true) const;

    /**
     * Check tree structure for parameter tree.
     * @param tree Parameter tree to check validity of.
     * @return Validty of tree. */
    virtual bool ValidParameterTree(TTree* tree) const;

    /**
     * Close MCMC output file. */
    void MCMCCloseOutputFile();

    /**
     * Marginalize from TTree.
     * @param autorange Flag for automatically choosing variable ranges to. */
    virtual void Remarginalize(bool autorange = true);

    /**
     * Update cholesky-decompositions for multivariate proposal function. */
    virtual bool UpdateCholeskyDecompositions();

    /** @} */

private:
    /**
     * Keeps variables that need to be both members and thread local.
     */
    struct MCMCThreadLocalStorage {
        /**
         * Store local proposal point  */
        std::vector<double> xLocal;

        /**
         * Random number generator */
        TRandom3* rng;

        /**
         * Constructor
         * @param dim Dimension of a temporary parameter vector */
        MCMCThreadLocalStorage(const unsigned& dim);

        /**
         * Copy constructor. */
        MCMCThreadLocalStorage(const MCMCThreadLocalStorage& other);

        /**
         * Assignment operator */
        MCMCThreadLocalStorage& operator = (const MCMCThreadLocalStorage& other);

        /**
         * Destructor */
        virtual ~MCMCThreadLocalStorage();
    };

    /**
     * Keep thread local variables private. */
    std::vector<MCMCThreadLocalStorage> fMCMCThreadLocalStorage;

    /**
     * Ensure that there are as many storages as chains */
    void SyncThreadStorage();

protected:

    /**
     * Set whether to fill 2D histogram y vs x: positive indices for
     * parameters; negative for observables, starting at -1 and going
     * more negative---observable index = -(index+1).
     * @param x Index of variable for horizontal axis.
     * @param y Index of variable for vertical axis.
     * @param flag Whether to fill 2D histogram. */
    void SetFillHistogram(int x, int y, bool flag);

    /**
     * Print model summary to stream.
     * @param ofi Output stream to print to. */
    virtual void PrintSummaryToStream(std::ofstream& ofi) const;

    /**
     * Print best fit to stream.
     * @param ofi Output stream to print to. */
    virtual void PrintBestFitToStream(std::ofstream& ofi) const;

    /**
     * Print marginalization to stream.
     * @param ofi Output stream to print to. */
    virtual void PrintMarginalizationToStream(std::ofstream& ofi) const;

    /**
     * Update Paramater TTree with scales and efficiencies. */
    void UpdateParameterTree();

    /**
     * return appropriate update interval
     * @param N total number of iterations to be accomplished.
     * @return Appropriate interval to output after. */
    unsigned  UpdateFrequency(unsigned N) const;

    /**
     * Name of the engine. */
    std::string fName;

    /**
     * Safe name of the engine for use in naming ROOT objects. */
    std::string fSafeName;

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
     * Number of iterations between scale adjustments and convergence
     * checks in pre-run. */
    unsigned fMCMCNIterationsPreRunCheck;

    /**
     * Number of iterations between clearing of convergence stats in pre-run.*/
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
    TFile* fMCMCOutputFile;

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
     * Minimum number of multivariate-proposal-function covariance updates to perform. */
    unsigned fMultivariateProposalFunctionUpdatesMinimum;

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

    /**
     * The phase of the run. */
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
     * flag for correcting R value for initial sampling variability. */
    bool fCorrectRValueForSamplingVariability;

    /**
     * The R-value criterion for convergence of parameters */
    double fMCMCRValueParametersCriterion;

    /** The R-values for each parameter */
    std::vector<double> fMCMCRValueParameters;

    /**
     * Random number generator */
    TRandom3* fRandom;

    /**
     * Vector of 1D marginalized distributions */
    std::vector<TH1*> fH1Marginalized;

    /**
     * Vector of 2D marginalized distributions. */
    std::vector<std::vector<TH2*> > fH2Marginalized;

    /**
     * Vector of pairs of indices for which 2D histograms should be
     * stored. a negative index indicates an observable, with
     * observable zero as -1, observable one as -2, etc. */
    std::vector<std::pair<int, int> > fRequestedH2;

    /**
     * The tree containing the Markov chains.*/
    TTree* fMCMCTree;

    /**
     * flag for whether MCMC Tree successfully loaded.*/
    bool fMCMCTreeLoaded;

    /**
     * flag for whether to reuse MCMC Tree's observables. */
    bool fMCMCTreeReuseObservables;

    unsigned int fMCMCTree_Chain;			///< (For writing to fMCMCTree) chain number
    unsigned int fMCMCTree_Iteration;	///< (For writing to fMCMCTree) iteration number
    double fMCMCTree_Prob;						///< (For writing to fMCMCTree) log(posterior)
    double fMCMCTree_LogLikelihood;		///< (For writing to fMCMCTree) log(likelihood)
    double fMCMCTree_LogPrior;				///< (For writing to fMCMCTree) log(prior)
    std::vector<double> fMCMCTree_Parameters;	///< (For writing to fMCMCTree) parameter set
    std::vector<double> fMCMCTree_Observables; ///< (For writing to fMCMCTree) observable set

    /**
     * The tree containing the parameter information.*/
    TTree* fParameterTree;

    /**
     * Vector of local modes. */
    std::vector<double> fLocalModes;

    /**
     * A BCH1D (with no histogram) for storing BCH1D drawing options. */
    BCH1D* fBCH1DdrawingOptions;

    /**
     * A BCH2D (with no histogram) for storing BCH2D drawing options. */
    BCH2D* fBCH2DdrawingOptions;

    /**
     * flag for rescaling of histograms after pre-run. */
    bool fRescaleHistogramRangesAfterPreRun;

    /**
     * factor for enlarging range of histograms when rescaling. */
    double fHistogramRescalePadding;

};

// ---------------------------------------------------------

#endif
