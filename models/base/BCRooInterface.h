#ifndef __BCROOINTERFACE__H
#define __BCROOINTERFACE__H

#include "../../BAT/BCModel.h"

#include <RooStats/MarkovChain.h>
#include <RooRealVar.h>

#include <list>
#include <utility>

class RooAbsReal;
class RooAbsData;
class RooAbsPdf;
class RooArgSet;
class RooArgList;

/** ---------------------------------------------------------
 * interface allowing to run BAT on a problem/data defined in a
 * standard RooFit workspace format
 * --------------------------------------------------------- */

class BCRooInterface : public BCModel
{
   public:

      BCRooInterface( const char* name = "", bool fillChain = false );

      ~BCRooInterface();

      /// Overloaded methods
      void DefineParameters();
      double LogAPrioriProbability(const std::vector<double> & parameters);
      double LogLikelihood(const std::vector<double> &parameters);

      /// Other method of this class
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
            const char* priorNuisanceName= "priorNuisance",
            const char* paramsName = "parameters",
            const char* listPOIName = "POI" );

      /// set the number of histogram bins for a specific parameter
      void SetNumBins(const char * parname, int nbins);
      /// set the number of histogram bins for all parameters
      void SetNumBins(int nbins);
      /// setup RooStats Markov Chain
      void SetupRooStatsMarkovChain();
      /// overloaded function from BCIntegrate to fill RooStats Markov Chain with every accepted step
      void MCMCIterationInterface();
      /// return the RooStats Markov Chain (empty if corresponding constructor option not set)
      RooStats::MarkovChain * GetRooStatsMarkovChain(){ return _roostatsMarkovChain;}
      RooArgSet* GetArgSetForMarkovChain(){return &_parametersForMarkovChainCurrent;}
      //RooArgSet GetArgSetForMarkovChainTest(){return _parametersForMarkovChain_test;}

   private:

      void AddToCurrentChainElement(double xij, int chainNum, int poiNum); //help function for construction of RooStats Markov Chain
      bool EqualsLastChainElement(int chainNum); //help function for construction of RooStats Markov Chain
      double GetWeightForChain(int chainNum); //help function for construction of RooStats Markov Chain

      RooAbsData* fData;        /// data to test
      RooAbsPdf*  fModel;       /// likelihood model describing the observables
      RooAbsReal* fNll;         /// pointer to negative log-likelihood function
      RooArgSet*  fObservables; /// list of observables measured for each event
      RooArgList* fParams;      /// list of parameters
      RooArgList* fParamsPOI;   /// list of parameters of interest
      RooAbsPdf*  fPrior;       /// function describing the prior probability of the parameters
      int _default_nbins;

      RooRealVar* priorhelpvar;
      bool _addeddummyprior;

      bool _fillChain;          ///decides if a RooStats Markov Chain is constructed along the way
      bool fFirstComparison;    /// checks if it is the first chain step (bookkeeping required for construction of the RooStats MarkovChain object)
      RooStats::MarkovChain * _roostatsMarkovChain;
      RooArgSet _parametersForMarkovChainPrevious; /// parameters of interest in previous step (bookkeeping required for construction of the RooStats MarkovChain object)
      RooArgSet _parametersForMarkovChainCurrent; /// parameters of interest in previous step (bookkeeping required for construction of the RooStats MarkovChain object)

      std::vector< std::vector<double> > fPreviousStep; /// bookkeeping required for construction of the RooStats MarkovChain object
      std::vector< std::vector<double> > fCurrentStep; /// bookkeeping required for construction of the RooStats MarkovChain object
      std::vector< double > fVecWeights; /// keep track of weights if proposal step not accepted (bookkeeping required for construction of the RooStats MarkovChain object)

      std::list< std::pair<const char*,int> > _nbins_list;
};

#endif
