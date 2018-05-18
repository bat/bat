#ifndef __BCROOINTERFACE__H
#define __BCROOINTERFACE__H

#include "../../BAT/BCModel.h"

#include <RooStats/MarkovChain.h>
#include <RooRealVar.h>

#include <list>
#include <string>
#include <utility>

class RooAbsReal;
class RooAbsData;
class RooAbsPdf;
class RooArgSet;
class RooArgList;

/**
 * Interface allowing to run BAT on a problem/data defined in a
 * standard RooFit workspace format.
 *
 * @warning The interface is available for backward compatibility but
 * it is essentially unmaintained. Known limitations are
 * + thread safety: this class is not thread safe.
 * + root 6 support: you cannot use BCRooInterface or BATCalculator from within cling. But you can use these two classes in compiled code.
 * + copy and assignment are untested
 */
class BCRooInterface : public BCModel
{
public:

    //! @nowarn
    BCRooInterface(const std::string& name = "", bool fillChain = false );

    ~BCRooInterface();
    //! @endnowarn

    /** \name Overloaded methods  */
    /** @{ */

    void DefineParameters();
    double LogAPrioriProbability(const std::vector<double>& parameters);
    double LogLikelihood(const std::vector<double>& parameters);
    /** @} */

    //! @nowarn
    void Initialize( RooAbsData& data,
                     RooAbsPdf& model,
                     RooAbsPdf& prior,
                     const RooArgSet* params,
                     const RooArgSet& listPOI );

    void Initialize( const char* rootFile,
                     const char* wsName = "batWS",
                     const char* dataName = "data",
                     const char* modelName = "model",
                     const char* priorName = "priorPOI",
                     const char* priorNuisanceName = "priorNuisance",
                     const char* paramsName = "parameters",
                     const char* listPOIName = "POI" );
    //! @endnowarn

    /** set the number of histogram bins for a specific parameter */
    void SetNumBins(const char* parname, int nbins);
    /** set the number of histogram bins for all parameters */
    void SetNumBins(int nbins);
    /** setup RooStats Markov Chain */
    void SetupRooStatsMarkovChain();
    /** overloaded function from BCIntegrate to fill RooStats Markov Chain with every accepted step */
    void MCMCIterationInterface();
    /** return the RooStats Markov Chain (empty if corresponding constructor option not set) */
    RooStats::MarkovChain* GetRooStatsMarkovChain() { return _roostatsMarkovChain;}
    //! @nowarn
    RooArgSet* GetArgSetForMarkovChain() {return &_parametersForMarkovChainCurrent;}
    //! @endnowarn
private:

    void AddToCurrentChainElement(double xij, int chainNum, int poiNum); ///< help function for construction of RooStats Markov Chain
    bool EqualsLastChainElement(int chainNum); ///< help function for construction of RooStats Markov Chain
    double GetWeightForChain(int chainNum); ///< help function for construction of RooStats Markov Chain

    RooAbsData* fData;        ///< data to test
    RooAbsPdf* fModel;        ///< likelihood model describing the observables
    RooAbsReal* fNll;         ///< pointer to negative log-likelihood function
    RooArgList* fParams;      ///< list of parameters
    RooArgList* fParamsPOI;   ///< list of parameters of interest
    RooAbsPdf* fPrior;        ///< function describing the prior probability of the parameters
    int _default_nbins;

    RooRealVar* priorhelpvar;
    bool _addeddummyprior;

    bool _fillChain;          ///<decides if a RooStats Markov Chain is constructed along the way
    bool fFirstComparison;    ///< checks if it is the first chain step (bookkeeping required for construction of the RooStats MarkovChain object)
    RooStats::MarkovChain* _roostatsMarkovChain;
    RooArgSet _parametersForMarkovChainPrevious; ///< parameters of interest in previous step (bookkeeping required for construction of the RooStats MarkovChain object)
    RooArgSet _parametersForMarkovChainCurrent; ///< parameters of interest in previous step (bookkeeping required for construction of the RooStats MarkovChain object)

    std::vector< std::vector<double> > fPreviousStep; ///< bookkeeping required for construction of the RooStats MarkovChain object
    std::vector< std::vector<double> > fCurrentStep; ///< bookkeeping required for construction of the RooStats MarkovChain object
    std::vector< double > fVecWeights; ///< keep track of weights if proposal step not accepted (bookkeeping required for construction of the RooStats MarkovChain object)

    std::list< std::pair<const char*, int> > _nbins_list;
};

#endif
