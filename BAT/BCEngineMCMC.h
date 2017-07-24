#ifndef __BCENGINEMCMC__H
#define __BCENGINEMCMC__H

/**
 * @class BCEngineMCMC
 * @brief An engine class for Markov Chain Monte Carlo
 * @author Daniel Kollar
 * @author Kevin Kr&ouml;ninger
 * @version 1.0
 * @date 08.2008
 * @details This class represents an engine class for Markov Chain
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

#include "BCH1D.h"
#include "BCH2D.h"
#include "BCLog.h"
#include "BCObservable.h"
#include "BCObservableSet.h"
#include "BCParameter.h"
#include "BCParameterSet.h"

#include <TMatrixD.h>
#include <TMatrixDSym.h>
#include <TRandom3.h>
#include <TVectorD.h>

#include <algorithm>
#include <limits>
#include <map>
#include <string>
#include <utility>
#include <vector>

// ROOT classes
class TF1;
class TFile;
class TH1;
class TH2;
class TTree;

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
    enum Phase {
        kPreRun     = -1, ///< In pre-run
        kUnsetPhase =  0, ///< Unset
        kMainRun    = +1  ///< In main run
    };

    /** An enumerator for markov-chain position initialization. */
    enum InitialPositionScheme {
        kInitContinue         = -1, ///< continue markov chain from last value (not yet available)
        kInitCenter           =  0, ///< select centers of parameter ranges
        kInitRandomUniform    =  1, ///< randomly distribute uniformly over parameter ranges
        kInitUserDefined      =  2, ///< initialize to user-provided points
        kInitRandomPrior      =  3  ///< randomly distribute according to factorized priors
    };

    /** @} */
    /** \name Structs */
    /** @{ */

    /** A struct for holding statistical information about samples. */
    struct Statistics {

        /** Constructor.
         * @param n_par number of parameters to calculate statistics for.
         * @param n_obs number of observables to calculate statistics for (sans efficiencies). */
        Statistics(unsigned n_par = 0, unsigned n_obs = 0);

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
        void Clear(bool clear_mode = true, bool clear_efficiency = true);

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

        /** addition assignment operator. */
        Statistics& operator += (const Statistics& rhs);

        /** update statistics given new values.
         * @param par Current parameter values.
         * @param obs Current user-define observables values. */
        void Update(double prob, const std::vector<double>& par, const std::vector<double>& obs);
    };

    /** @} */
    /** \name Constructors and destructor */
    /** @{ */

    /**
     * Default constructor. */
    BCEngineMCMC(const std::string& name = "model");

    /**
     * Copy constructor. */
    BCEngineMCMC(const BCEngineMCMC& enginemcmc);

    /**
     * Read in MCMC constructor.
     * @param filename Path of file holding model.
     * @param name Name of model (file should contain TTree's [name]_mcmc and [name]_parameters.\n
     * if empty string is given, properly matching TTrees are searched for in the file.
     * @param loadObservables Flag for whether to load observables for file (true; default) or to let user respecify observables.*/
    BCEngineMCMC(const std::string& filename, const std::string& name, bool loadObservables = true);

    /**
     * Destructor. */
    virtual ~BCEngineMCMC();

    /** @} */
    /** \name swap */
    /** @{ */

    /** swap */
    friend void swap(BCEngineMCMC& A, BCEngineMCMC& B);

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
    unsigned GetNChains() const
    { return fMCMCNChains; }

    /**
     * @return lag of the Markov chains */
    unsigned GetNLag() const
    { return fMCMCNLag; }

    /**
     * @return number of iterations */
    const std::vector<unsigned>& GetNIterations() const
    { return fMCMCNIterations; }

    /**
     * @return current iterations */
    int GetCurrentIteration() const
    { return fMCMCCurrentIteration; }

    /**
     * @return current chain index. If no chain is currently selected
     * (for example, not within mcmc), return 0. */
    unsigned GetCurrentChain() const;

    /**
     * @return number of iterations needed for all chains to
     * converge simultaneously. A value of -1 indicates that the chains did not converge. */
    int GetNIterationsConvergenceGlobal() const
    { return fMCMCNIterationsConvergenceGlobal; }

    /**
     * @return the number of iterations the prerun actually took. If
     * the prerun wasn't run, that's 0. Else it is bounded by
     * GetNIterationsPreRunMin() and GetNIterationsPreRunMax().
     */
    unsigned GetNIterationsPreRun() const;

    /**
     * @return minimum number of pre-run iterations for a Markov chain */
    unsigned GetNIterationsPreRunMin() const
    { return fMCMCNIterationsPreRunMin; }

    /**
     * @return maximum number of pre-run iterations for a Markov chain */
    unsigned GetNIterationsPreRunMax() const
    { return fMCMCNIterationsPreRunMax; }

    /**
     * @return number of iterations for a Markov chain */
    unsigned GetNIterationsRun() const
    { return fMCMCNIterationsRun; }

    /**
     * @return number of iterations between scale adjustments and convergence checking during pre-run. */
    unsigned GetNIterationsPreRunCheck() const
    { return fMCMCNIterationsPreRunCheck; }

    /**
     * @return number of iterations between clearings of statistics for convergence checking. */
    unsigned GetPreRunCheckClear()
    { return fMCMCPreRunCheckClear; }

    /**
     * @return minimum efficiency required for a chain. */
    double GetMinimumEfficiency() const
    { return fMCMCEfficiencyMin; }

    /**
     * @return maximum efficiency required for a chain. */
    double GetMaximumEfficiency() const
    { return fMCMCEfficiencyMax; }

    /**
     * @return proposal function scale factor lower limit */
    double GetScaleFactorLowerLimit() const
    { return fMCMCScaleFactorLowerLimit; }

    /**
     * @return proposal function scale factor upper limit */
    double GetScaleFactorUpperLimit() const
    { return fMCMCScaleFactorUpperLimit; }

    /**
     * @return scale factor for all parameters and chains */
    const std::vector<std::vector<double> >& GetScaleFactors() const
    { return fMCMCProposalFunctionScaleFactor; }

    /**
     * @return current point of each Markov chain */
    const std::vector<std::vector<double> >& Getx() const
    { return fMCMCx; }

    /**
     * @param c index of the Markov chain
     * @return current point of the Markov chain */
    const std::vector<double>& Getx(unsigned c) const
    { return fMCMCx.at(c); }

    /**
     * @param c chain index
     * @param p parameter index
     * @return parameter of the Markov chain */
    double Getx(unsigned c, unsigned p) const
    { return fMCMCx.at(c).at(p); }

    /**
     * @return log of the probability of the current points of each Markov chain */
    const std::vector<double>& GetLogProbx() const
    { return fMCMCprob; }

    /**
     * @return log of the probability of the current points of the Markov chain.
     * @param c chain index */
    double GetLogProbx(unsigned c) const
    { return fMCMCprob.at(c); }

    /**
     * @return pointer to the phase of a run. */
    BCEngineMCMC::Phase GetPhase() const
    { return fMCMCPhase; }

    /**
     * @return flag which defined initial position */
    BCEngineMCMC::InitialPositionScheme GetInitialPositionScheme() const
    { return fInitialPositionScheme; }

    /**
     * @return maximum number of allowed attempts to set the initial position */
    unsigned GetInitialPositionAttemptLimit() const
    { return fInitialPositionAttemptLimit; }

    /**
     * @return whether to use a multivariate proposal function. */
    bool GetProposeMultivariate() const
    { return fMCMCProposeMultivariate; }

    /**
     * @return Degree of freedom of multivariate proposal
     * function. Anything >0 is the degree of freedom of a
     * multivariate Student's t distribution, <= 0 is a multivariate
     * Gaussian */
    double GetProposalFunctionDof() const
    { return fMCMCProposalFunctionDof; }

    /**
     * @return number of updates to multivariate-proposal-function covariance performed. */
    unsigned GetMultivariateCovarianceUpdates() const
    { return fMultivariateCovarianceUpdates; }

    /*
     * @return weighting parameter for multivariate proposal function covariance update. */
    double GetMultivariateCovarianceUpdateLambda() const
    { return fMultivariateCovarianceUpdateLambda; }

    /**
     * @return multivariate-proposal-function Cholesky decomposition nudge size. */
    double GetMultivariateEpsilon() const
    { return fMultivariateEpsilon; }

    /**
     * @return multivariate-proposal-function scale multiplier. */
    double GetMultivariateScaleMultiplier() const
    { return fMultivariateScaleMultiplier; }

    /**
     * @return R-value criterion for parameters */
    double GetRValueParametersCriterion() const
    { return fMCMCRValueParametersCriterion; }

    /**
     * @return vector of R values for parameters */
    const std::vector<double>& GetRValueParameters() const
    { return fMCMCRValueParameters; }

    /**
     * @return R-value for a parameter
     * @param index parameter index */
    double GetRValueParameters(unsigned index) const
    { return fMCMCRValueParameters.at(index); }

    /** Flag for correcting convergence checking for initial sampling variability. */
    bool GetCorrectRValueForSamplingVariability() const
    { return fCorrectRValueForSamplingVariability; }

    /**
     * @return the flag if MCMC has been performed or not */
    bool GetFlagRun() const
    { return fMCMCFlagRun; }

    /**
     * @return the flag if MCMC pre-run should be performed or not */
    bool GetFlagPreRun() const
    { return fMCMCFlagPreRun; }

    /**
     * Retrieve the tree containing the Markov chain. */
    TTree* GetMarkovChainTree() const
    { return fMCMCTree;}

    /**
     * Retrieve the tree containing the parameter information. */
    TTree* GetParameterTree() const
    { return fParameterTree;}

    /**
     * Retrieve output file for MCMC. */
    TFile* GetOutputFile() const
    { return fMCMCOutputFile; }

    /**
     * Get combined statistics for all chains. */
    const BCEngineMCMC::Statistics& GetStatistics() const
    { return fMCMCStatistics_AllChains; }

    /**
     * Get MCMC statistics for one chain.
     * @param c Chain to get statistics of. */
    const BCEngineMCMC::Statistics& GetStatistics(unsigned c) const
    { return fMCMCStatistics.at(c); }

    /**
     * Get vector of MCMC statistics for each chain separately. */
    const std::vector<BCEngineMCMC::Statistics>& GetStatisticsVector() const
    { return fMCMCStatistics; }

    /**
     * @return Flag for whether to rescale histogram ranges to fit MCMC reach after pre-run. */
    bool GetRescaleHistogramRangesAfterPreRun() const
    { return fRescaleHistogramRangesAfterPreRun; }

    /**
     * @return Factor by which to enlarge histogram ranges when rescaling to add padding beyond range. */
    double GetHistogramRescalePadding() const
    { return fHistogramRescalePadding; }

    /**
     * @return vector indices for order of printing 1D histograms. */
    virtual std::vector<unsigned> GetH1DPrintOrder() const;

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
    TH1* GetMarginalizedHistogram(const std::string& name) const
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
    TH2* GetMarginalizedHistogram(const std::string& name1, const std::string& name2) const
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
    BCH1D GetMarginalized(const std::string& name) const
    { return GetMarginalized(fParameters.Index(name)); }

    /**
     * Obtain the individual marginalized distributions
     * with respect to one parameter.
     * @note Ownership of the returned heap object is conferred to the caller.
     * @param index The parameter index
     * @return 1D marginalized probability */
    BCH1D GetMarginalized(unsigned index) const;

    /**
     * Obtain the individual marginalized distributions
     * with respect to two parameters.
     * @note The most efficient method is to access by indices.
     * @note Ownership of the returned heap object is conferred to the caller.
     * @param name1 Name of first parameter
     * @param name2 Name of second parameter
     * @return 2D marginalized probability */
    BCH2D GetMarginalized(const std::string& name1, const std::string& name2) const
    { return GetMarginalized(fParameters.Index(name1), fParameters.Index(name2)); }

    /**
     * Obtain the individual marginalized distributions
     * with respect to two parameters.
     * @note Ownership of the returned heap object is conferred to the caller.
     * @param index1 Index of first parameter
     * @param index2 Index of second parameter
     * @return 2D marginalized probability */
    BCH2D GetMarginalized(unsigned index1, unsigned index2) const;

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
     * @deprecated Instead call GetParameters().Get(name)
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
    virtual const std::vector<double>& GetBestFitParameters() const
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
    bool GetReuseObservables() const
    { return fMCMCTreeReuseObservables; }

    /**
     * @return BCH1D object that stores drawing options for all BCH1D's. */
    BCH1D& GetBCH1DdrawingOptions()
    { return fBCH1DdrawingOptions; }

    /**
     * @return BCH2D object that stores drawing options for all BCH2D's. */
    BCH2D& GetBCH2DdrawingOptions()
    { return fBCH2DdrawingOptions; }

    /** @} */
    /** \name Setters */
    /** @{ */

    /**
     * Sets the name of the engine.
     * @param name Name of the engine */
    void SetName(const std::string& name);

    /**
     * Set scale factor lower limit */
    void SetScaleFactorLowerLimit(double l)
    { fMCMCScaleFactorLowerLimit = l; }

    /**
     * Set scale factor upper limit */
    void SetScaleFactorUpperLimit(double l)
    { fMCMCScaleFactorUpperLimit = l; }

    /**
     * Set the initial scale factors for the factorized proposal function.
     * The same factors are used for every chain and get updated during the pre-run.
     * @param scale a vector of doubles containing the scale factors for each parameter. */
    void SetInitialScaleFactors(const std::vector<double>& scale)
    { fMCMCInitialScaleFactors = scale; }

    /**
     * Sets the number of Markov chains which are run in parallel. */
    void SetNChains(unsigned n)
    { fMCMCNChains = n; }

    /**
     * Sets the lag of the Markov chains */
    void SetNLag(unsigned n);

    /**
     * Sets the maximum number of iterations in the pre-run. */
    void SetNIterationsPreRunMax(unsigned n)
    { fMCMCNIterationsPreRunMax = n; }

    /**
     * Sets the number of iterations. */
    void SetNIterationsRun(unsigned n)
    { fMCMCNIterationsRun = n; }

    /**
     * Sets the minimum number of iterations in the pre-run */
    void SetNIterationsPreRunMin(unsigned n)
    { fMCMCNIterationsPreRunMin = n; }

    /**
     * Sets the number of iterations between scale adjustments and
     * convergence checks in the pre-run. */
    void SetNIterationsPreRunCheck(unsigned n)
    { fMCMCNIterationsPreRunCheck = n; }

    /**
     * Sets the number of prerun checks to make inbetween statistics clearing. */
    void SetPreRunCheckClear(unsigned n)
    { fMCMCPreRunCheckClear = n; }

    /**
     * Sets the minimum efficiency required for a chain. */
    void SetMinimumEfficiency(double efficiency)
    { fMCMCEfficiencyMin = efficiency; }

    /**
     * Sets the maximum efficiency required for a chain. */
    void SetMaximumEfficiency(double efficiency)
    { fMCMCEfficiencyMax = efficiency; }

    /**
     * Set the random number seed */
    void SetRandomSeed(unsigned seed);

    /**
     * Sets the initial positions for all chains.
     * @param x0s initial positions for all chains. */
    void SetInitialPositions(const std::vector<double>& x0s);

    /**
     * Sets the initial positions for all chains.
     * @param x0s initial positions for all chains. */
    void SetInitialPositions(const std::vector< std::vector<double> >& x0s)
    { fMCMCInitialPosition = x0s; SetInitialPositionScheme(BCEngineMCMC::kInitUserDefined); }

    /**
     * Sets flag which defines initial position.  */
    void SetInitialPositionScheme(BCEngineMCMC::InitialPositionScheme scheme)
    { fInitialPositionScheme = scheme; }

    /**
     * Sets maximum number of attempts to find a valid initial position. */
    void SetInitialPositionAttemptLimit(unsigned n)
    { fInitialPositionAttemptLimit = n; }

    /**
     * Set `flag` to `true` to turn on the multivariate proposal for
     * MCMC based on (Haario et al., 2001) where the covariance is
     * learned from the prerun.
     *
     * If `flag == false`, use a factorized proposal in which every
     * parameter is varied individually, one after the other. This
     * means for N parameters, N calls to the target density occur
     * until every parameter has been (attempted to be) varied exactly
     * once. In contrast, with the multivariate proposal all
     * parameters are varied simultaneously and a move can occur for a
     * single call to the target density.
     *
     * For both multivariate and factorized proposal, the acceptance
     * rate of the Markov chain is tuned to lie within the limits
     * given by SetMinimumEfficiency() and SetMaximumEfficiency().*/
    void SetProposeMultivariate(bool flag)
    {
        fMCMCProposeMultivariate = flag;
    }

    /**
     * Set the degree of freedom of the proposal function for MCMC.
     *
     * The default `dof == 1` represents a Cauchy distribution. For
     * any positive value of `dof`, a Student's t distribution with
     * the corresponding degree of freedom is used. For `dof <= 0`,
     * a Gaussian is used.
     *
     * A small positive degree of freedom leads to fat tails in the
     * proposal. This makes it easier to make a large jump in a single
     * iteration but generally leads to a lower acceptance rate being optimal.
     */
    void SetProposalFunctionDof(double dof = 1)
    {
        fMCMCProposalFunctionDof = dof;
    }

    /**
     * Set weighting for multivariate proposal function covariance update.
     * value forced into [0, 1] */
    void SetMultivariateCovarianceUpdateLambda(double l)
    { fMultivariateCovarianceUpdateLambda = std::max<double>(0, std::min<double>(1, l)); }

    /**
     * Sets multivariate-proposal-function cholesky-decomposition nudge. */
    void SetMultivariateEpsilon(double epsilon)
    { fMultivariateEpsilon = epsilon; }

    /**
     * Sets multivariate-proposal-function scale multiplier. */
    void SetMultivariateScaleMultiplier(double s)
    { fMultivariateScaleMultiplier = s; }

    /** Sets whether to fill histograms.
     * Applies to 1D and 2D histogams of all parameters and observables so far added. */
    void SetFlagFillHistograms(bool flag)
    { fParameters.FillHistograms(flag); fObservables.FillHistograms(flag); }

    /** Sets the whether to fill histograms.
     * Applies to all parameters and observables so far added. */
    void SetFlagFillHistograms(bool flag_1d, bool flag_2d)
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
    void SetFillHistogramParPar(const std::string& x, const std::string& y, bool flag = true)
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
    void SetFillHistogramParObs(const std::string& x, const std::string& y, bool flag = true)
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
    void SetFillHistogramObsObs(const std::string& x, const std::string& y, bool flag = true)
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
    void SetFillHistogramObsPar(const std::string& x, const std::string& y, bool flag = true)
    { SetFillHistogramObsPar(fObservables.Index(x), fParameters.Index(y), flag); }

    /** Sets the flag if a prerun should be performed or not. */
    void SetFlagPreRun(bool flag)
    { fMCMCFlagPreRun = flag; }

    /**
     * Sets the parameter R-value criterion for convergence of all chains */
    void SetRValueParametersCriterion(double r)
    { fMCMCRValueParametersCriterion = r; }

    /** Set flag to correct convergence checking for initial sampling variability. */
    void SetCorrectRValueForSamplingVariability(bool flag = true)
    { fCorrectRValueForSamplingVariability = flag; }

    /**
     * Set the precision for the MCMC run. */
    void SetPrecision(BCEngineMCMC::Precision precision);

    /**
     * Copy precision for the MCMC run from other model. */
    void SetPrecision(const BCEngineMCMC* other)
    { if (other) SetPrecision(*other); }

    /**
     * Copy precision for the MCMC run from other model. */
    void SetPrecision(const BCEngineMCMC& other);

    /**
     * Set the number of bins for the marginalized distribution of all parameters
     * @param nbins Number of bins */
    void SetNbins(unsigned int nbins)
    { fParameters.SetNBins(nbins); fObservables.SetNBins(nbins); }

    /**
     * @param flag Flag whether to reuse user-defined observables from MCMC tree when looping through it. */
    void SetReuseObservables(bool flag)
    { fMCMCTreeReuseObservables = flag; }

    /**
     * Set flag for rescaling histogram ranges after pre-run. */
    void SetRescaleHistogramRangesAfterPreRun(bool flag = true)
    { fRescaleHistogramRangesAfterPreRun = flag; }

    /**
     * Set enlargement factor of range for when rescaling. */
    void SetHistogramRescalingPadding(double factor)
    { fHistogramRescalePadding = factor; }

    /**
     * Turn on/off writing of Markov chain to root file.
     * If setting true, you must first set filename with function with filename arguments.
     * @param flag Flag for writing Markov chain (run and prerun) to ROOT file. */
    void WriteMarkovChain(bool flag)
    { WriteMarkovChainRun(flag); WriteMarkovChainPreRun(flag); }

    /**
     * Turn on/off writing of Markov chain to root file during run.
     * If setting either true, you must first set filename with function with filename arguments.
     * @param flag Flag for writing run Markov chain to ROOT file (true) or not (false). */
    void WriteMarkovChainRun(bool flag);

    /**
     * Turn on/off writing of Markov chain to root file during prerun.
     * If setting either true, you must first set filename with function with filename arguments.
     * @param flag Flag for writing prerun Markov chain to ROOT file (true) or not (false). */
    void WriteMarkovChainPreRun(bool flag);

    /**
     * Turn on writing of Markov chain to root file.
     * @param filename Name of file to write chain to.
     * @param option file-open options (TFile), must be "NEW", "CREATE", "RECREATE", or "UPDATE" (i.e. writeable).
     * @param flag_run Flag for writing run Markov chain to ROOT file (true) or not (false).
     * @param flag_prerun Flag for writing prerun Markov chain to ROOT file (true) or not (false). */
    void WriteMarkovChain(const std::string& filename, const std::string& option, bool flag_run = true, bool flag_prerun = true);

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
    void SetPriorConstant(const std::string& name)
    { SetPriorConstant(fParameters.Index(name)); }

    /**
     * @deprecated Instead call: GetPrarameter(index)->SetPrior(new BCTF1Prior(f))
     * Set prior for a parameter.
     * @param index The parameter index
     * @param f A function describing the prior
     * @param logL Whether function is of log of prior (true) or just prior (false)*/
    void SetPrior(unsigned index, TF1& f, bool logL = true);

    /**
     * @deprecated Instead call: GetParameter(name)->SetPrior(new BCTF1Prior(f))
     * Set prior for a parameter.
     * @param name The parameter name
     * @param f A function describing the prior
     * @param logL Whether function is of log of prior (true) or just prior (false)*/
    void SetPrior(const std::string& name, TF1& f, bool logL = true)
    { return SetPrior(fParameters.Index(name), f, logL); }

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
    void SetPriorDelta(const std::string& name, double value)
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
    void SetPriorGauss(const std::string& name, double mean, double sigma)
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
    void SetPriorGauss(const std::string& name, double mode, double sigma_below, double sigma_above)
    {	SetPriorGauss(fParameters.Index(name), mode, sigma_below, sigma_above); }

    /**
     * @deprecated Instead call: GetParameter(index)->SetPrior(new BCTH1Prior(h,interpolate))
     * Set prior for a parameter.
     * @param index parameter index
     * @param h histogram describing the prior
     * @param interpolate whether or not to use linear interpolation */
    void SetPrior(unsigned index, TH1& h, bool interpolate = false);

    /**
     * @deprecated Instead call: GetParameter(name)->SetPrior(new BCTH1Prior(h,interpolate))
     * Set prior for a parameter.
     * @param name parameter name
     * @param h histogram describing the prior
     * @param interpolate whether or not to use linear interpolation
     * @return success of action. */
    void SetPrior(const std::string& name, TH1& h, bool interpolate = false)
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
     * @param option Options passed to ROOT's TFile::Open.
     * @param closeExistingFile if file is already open, whether to close it after writing. */
    void WriteMarginalizedDistributions(const std::string& filename, const std::string& option, bool closeExistingFile = false);

    /**
     * Prints a summary to the logs. */
    virtual void PrintSummary() const;

    /**
     * Print parameters
     * @param P vector of the parameter values to be printed
     * @param output pointer to the output function to be used, which defaults to BCLog::OutSummary */
    void PrintParameters(const std::vector<double>& P, void (*output)(const std::string&) = BCLog::OutSummary) const;

    /**
     * Print all marginalizations.
     * @param filename Path to file to print to
     * @param hdiv Number of columns of plots per page
     * @param vdiv Number of rows of plots per page
     * @return Number of plots printed */
    unsigned PrintAllMarginalized(const std::string& filename, unsigned hdiv = 1, unsigned vdiv = 1) const;

    /**
     * Print a summary plot for the parameters and user-defined observables.
     * @param filename Path to filename to print to.
     * @param npar Number of parameters per page, print all on one page if set to zero or negative
     * @param interval_content Probability mass to display in smallest X interval band
     * @param quantile_values Vector of quantile values to draw
     * @param rescale_ranges Flag for rescaling to range surveyed by MCMC chains
     * @return Number of pages printed. */
    unsigned PrintParameterPlot(const std::string& filename, int npar = 10, double interval_content = 68e-2, std::vector<double> quantile_values = std::vector<double>(0), bool rescale_ranges = true) const;

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
    bool PrintCorrelationMatrix(const std::string& filename = "matrix.pdf") const;

    /**
     * Print a correlation plot for the parameters.
     * @param filename Path to file to print correlation plot to.
     * @param include_observables Flag for including observables (default: true)
     * @return Success of action. */
    bool PrintCorrelationPlot(const std::string& filename = "correlation.pdf", bool include_observables = true) const;

    /**
     * Print a LaTeX table of the parameters.
     * @param filename Path to file tp print LaTeX table of parameters to.
     * @return Success of action. */
    bool PrintParameterLatex(const std::string& filename) const;

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
     * @deprecated Instead use GetParameters().Add(...)
     * Adds a parameter.
     * @param name Name of parameter
     * @param min minimum value of the parameter
     * @param max maximum value of the parameter
     * @param latexname Optional latexname used for plotting
     * @param unitstring Unit string to be printed for parameter.
     * @return Success of action. */
    virtual bool AddParameter(const std::string& name, double min, double max, const std::string& latexname = "", const std::string& unitstring = "")
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
    virtual bool AddObservable(const std::string& name, double min, double max, const std::string& latexname = "", const std::string& unitstring = "")
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
     * to calculate user observables.
     * @param pars Parameter set to evaluate observables for. */
    virtual void CalculateObservables(const std::vector<double>& /*pars*/)
    {}

    /**
     * The default proposal function is a
     * Breit-Wigner random walk. It can be overloaded by the user to set the proposal
     * function.
     * @param ichain the chain index
     * @param ipar the parameter index
     * @return the unscaled proposal point */
    virtual double ProposalFunction(unsigned ichain, unsigned ipar);

    /**
     * Returns a proposal point for the Metropolis algorithm.
     * @param chain chain index
     * @param x proposal point
     * @return flag indicating whether the new point lies within the allowed range */
    bool GetProposalPointMetropolis(unsigned chain, std::vector<double>& x);

    /**
     * Returns a proposal point for the Metropolis algorithm.
     * @param chain chain index
     * @param x proposal point
     * @return flag indicating whether the new point lies within the allowed range */
    bool GetProposalPointMetropolis(unsigned chain, unsigned parameter, std::vector<double>& x);

    /**
     * Generates a new point using the Metropolis algorithm for all chains.
     * @return Whether proposed point was accepted (true) or previous point was kept (false). */
    bool GetNewPointMetropolis();

    /** Generates a new point using the Metropolis algorithm for one chain.
     * @param chain chain index
     * @return Whether proposed point was accepted (true) or previous point was kept (false). */
    bool GetNewPointMetropolis(unsigned chain);

    /** Generates a new point using the Metropolis algorithm for one chain, varying only one parameter's value
     * @param chain chain index
     * @param parameter index of parameter to vary
     * @return Whether proposed point was accepted (true) or previous point was kept (false). */
    bool GetNewPointMetropolis(unsigned chain, unsigned parameter);

    /**
     * Updates statistics: fill marginalized distributions */
    void InChainFillHistograms();

    /**
     * Updates statistics: write chains to file */
    void InChainFillTree();

    /**
     * Needs to be overloaded in the derived class.
     * @param parameters Parameter set to evaluate at.
     * @return natural logarithm of the function to map with MCMC */
    virtual double LogEval(const std::vector<double>& parameters) = 0;

    /**
     * Runs Metropolis algorithm.
     * @return Success of action. */
    bool Metropolis();

    /**
     * Runs a pre run for the Metropolis algorithm.
     * @return Success of action. */
    bool MetropolisPreRun();

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
    static double RValue(const std::vector<double>& means, const std::vector<double>& variances, unsigned n, bool correctForSamplingVariability = true);

    /**
     * Resets all containers used in MCMC and initializes starting points. */
    void MCMCInitialize();

    /**
     * User hook called from MCMCInitialize().
     *
     * MCMCUserInitialize() is called after all settings for the
     * upcoming MCMC run are fixed but before the initial point is
     * chosen and before any call to user methods such as
     * LogLikelihood() or LogAPrioriProbability() are
     * issued. MCMCUserInitialize() is useful for example to allocate
     * separate copies of objects within the user model, one per
     * chain, for thread safety.
     *
     * @note Any error inside MCMCUserInitialize() should be signaled via an exception. */
    virtual void MCMCUserInitialize()
    {}

    /**
     * Reset the MCMC variables. */
    virtual void ResetResults();

    /**
     * Interface allowing to execute arbitrary code for each iteration
     * of the MCMC while running the chains after applying the lag but
     * before filling histograms or the output tree.
     *
     * With Remarginalize(), there is no lag, and MCMCUserIterationInterface()
     * is called on every iteration before filling histograms.
     *
     * In both instances of use, CalculateObservables() has already been called.
     *
     * This function is not necessarily called for each call to the likelihood.
     * Its frequency is affected by the lag; and for the factorized proposal function,
     * it is called only after _all_ parameters have been updated.
     *
     * This method needs to be overloaded by the user to do anything.
     *
     * @note For uncertainty propagation, it is much more convenient BCObservable and CalculateObservables() */
    virtual void MCMCUserIterationInterface()
    {}

    /**
     * Interface allowing to execute arbitrary code for each new point
     * of the MCMC whether it is accepted or not. This method needs to be overloaded in the derived
     * class
     *
     * @note This method is called for every call to the likelihood.
     *
     * @param point point that was generated and checked
     * @param ichain index of the chain
     * @param accepted flag whether or not the point was accepted for the chain */
    virtual void MCMCCurrentPointInterface(const std::vector<double>& /*point*/, int /*ichain*/, bool /*accepted*/)
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
     * @return Whether tree's parameters match those of model. */
    virtual bool ParameterTreeMatchesModel(TTree* partree, bool checkObservables = true);

    /**
     * Load previous MCMC run.
     * @param filename Pathname of file containing model trees.
     * @param mcmcTreeName Name of tree inside file containing MCMC, empty string (default) loads [modelname]_mcmc.
     * @param parameterTreeName Name of tree inside file containing parameter list, empty string (default) loads [modelname]_parameters.
     * @param loadObservables Flag for whether to load observables from parameter list and MCMC trees. */
    virtual bool LoadMCMC(const std::string& filename, const std::string& mcmcTreeName = "", const std::string& parameterTreeName = "", bool loadObservables = true);

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
     * Close the root output file. */
    void CloseOutputFile();

    /**
     * Marginalize from TTree.
     * @param autorange Flag for automatically choosing variable ranges to. */
    virtual void Remarginalize(bool autorange = true);

    /**
     * Update cholesky-decompositions for multivariate proposal function. */
    virtual bool UpdateCholeskyDecompositions();

    /**
     * Keep track of which chain is currently computed (within a thread).
     *
     * @warning Call this method only if you know what you are doing!
     */
    void UpdateChainIndex(int chain);

    /** @} */

private:
    /**
     * Keeps variables that need to be both members and thread local.
     */
    struct ThreadLocalStorage {
        /**
         * Store local proposal point  */
        std::vector<double> xLocal;

        /**
         * Random number generator */
        TRandom3* rng;

        /**
         * Temp vector for matrix multiplication in multivariate proposal */
        TVectorD yLocal;

        /**
         * Constructor
         * @param dim Dimension of a temporary parameter vector
         */
        ThreadLocalStorage(unsigned dim);

        /**
         * Copy constructor. */
        ThreadLocalStorage(const ThreadLocalStorage& other);

        /**
         * Assignment operator */
        ThreadLocalStorage& operator = (ThreadLocalStorage other);

        /** swap, can't be a friend this time because this struct is private */
        void swap(BCEngineMCMC::ThreadLocalStorage& A, BCEngineMCMC::ThreadLocalStorage& B);

        /**
         * Destructor */
        virtual ~ThreadLocalStorage();

        /**
         * Scale a Gaussian random variate such that it becomes a student's t variate;
         * i.e. return `sqrt(dof / chi2)`, where `chi2` is a variate from a chi2 distribution
         * with `dof` degrees of freedom.
         *
         * If `dof <= 0`, simply return 1 to keep the Gaussian distribution.
         */
        double scale(double dof);
    };

    /**
     * Keep thread local variables private. */
    std::vector<ThreadLocalStorage> fMCMCThreadLocalStorage;

    /**
     * Ensure that there are as many storages as chains */
    void SyncThreadStorage();

    typedef std::map<int, unsigned> ChainIndex_t;
    ChainIndex_t fChainIndex;

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
     * Print model summary to log. */
    virtual void PrintModelSummary() const;

    /**
     * Print best fit to log. */
    virtual void PrintBestFitSummary() const;

    /**
     * Get string summarizing best fit for single variable.
     * @param i Index of variable to summarize. */
    virtual std::string GetBestFitSummary(unsigned i) const;

    /**
     * Print marginalization to log. */
    virtual void PrintMarginalizationSummary() const;

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
     * Number of iterations between scale adjustments and convergence
     * checks in pre-run. */
    unsigned fMCMCNIterationsPreRunCheck;

    /**
     * Number of iterations between clearing of convergence stats in pre-run.*/
    unsigned fMCMCPreRunCheckClear;

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

    /**
     * Output file for for writing MCMC Tree. */
    TFile* fMCMCOutputFile;

    /**
     * Output filename for for writing MCMC Tree. */
    std::string fMCMCOutputFilename;

    /**
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
    std::vector<std::vector<double> > fMCMCProposalFunctionScaleFactor;

    /**
     * User-provided initial values of the scale factors of the factorized proposal function. */
    std::vector<double> fMCMCInitialScaleFactors;

    /**
     * Covariance matrices used in multivariate proposal functions. */
    std::vector<TMatrixDSym> fMultivariateProposalFunctionCovariance;

    /**
     * Cholesky decompositions for multivariate proposal function.
     * Index is over chain. */
    std::vector<TMatrixD> fMultivariateProposalFunctionCholeskyDecomposition;

    /**
     * Number of multivariate-proposal-function covariance updates performed. */
    unsigned fMultivariateCovarianceUpdates;

    /**
     * weighting parameter for multivariate-proposal-function covariance update. */
    double fMultivariateCovarianceUpdateLambda;

    /**
     * multivariate-proposal-function cholesky-decomposition nudge. */
    double fMultivariateEpsilon;

    /**
     * factor to multiply or divide scale factors by in adjusting multivariate-proposal-function scales. */
    double fMultivariateScaleMultiplier;

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
     * The minimum required efficiency for MCMC */
    double fMCMCEfficiencyMin;

    /**
     * The maximum allowed efficiency for MCMC */
    double fMCMCEfficiencyMax;

    /**
     * Variable which defines the initial position.
     * See enum MCMCInitialPosition for possible values. */
    BCEngineMCMC::InitialPositionScheme fInitialPositionScheme;

    /**
     * Maximum number of attempts to make to set the initial position. */
    unsigned fInitialPositionAttemptLimit;

    /**
     * Flag for using multivariate proposal function. */
    bool fMCMCProposeMultivariate;

    /**
     * Degree of freedom of Student's t proposal. If <= 0, use Gaussian proposal.
     */
    double fMCMCProposalFunctionDof;

    /**
     * The phase of the run. */
    BCEngineMCMC::Phase fMCMCPhase;

    /**
     * The current points of each Markov chain. */
    std::vector<std::vector<double> > fMCMCx;

    /**
     * The current values of the user-defined observables for each
     * Markov chain. */
    std::vector<std::vector<double> > fMCMCObservables;

    /**
     * Statistics for each Markov chain. */
    std::vector<BCEngineMCMC::Statistics> fMCMCStatistics;

    /**
     * Statistics across all Markov chains. */
    BCEngineMCMC::Statistics fMCMCStatistics_AllChains;

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
    TRandom3 fRandom;

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
    BCH1D fBCH1DdrawingOptions;

    /**
     * A BCH2D (with no histogram) for storing BCH2D drawing options. */
    BCH2D fBCH2DdrawingOptions;

    /**
     * flag for rescaling of histograms after pre-run. */
    bool fRescaleHistogramRangesAfterPreRun;

    /**
     * factor for enlarging range of histograms when rescaling. */
    double fHistogramRescalePadding;

};

// ---------------------------------------------------------

#endif
