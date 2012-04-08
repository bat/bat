#include <RooGlobalFunc.h>
#include <RooMsgService.h>
#include <RooProdPdf.h>
#include <RooWorkspace.h>
#include <TFile.h>

#include "../../BAT/BCMath.h"
#include <iostream>

#include "BCRooInterface.h"

#include <RooUniform.h>

//for testing
#include "RooRealVar.h"
#include "RooAbsReal.h"
#include "RooRandom.h"




// ---------------------------------------------------------
void BCRooInterface::Initialize( RooAbsData& data,
             RooAbsPdf& model,
             RooAbsPdf& prior_trans,
             const RooArgSet* params,
             const RooArgSet& listPOI )
{

  fData = &data;
  fModel = &model;

  // make the product of both priors to get the full prior probability function
  RooAbsPdf* prior_total = &prior_trans;
  if (prior_total!=0 ) {
    fPrior = prior_total;
  }
  else {
         std::cout << "No prior PDF: without taking action the program would crash\n";
         std::cout << "No prior PDF: adding dummy uniform pdf on the interval [0..1]\n";
         priorhelpvar = new RooRealVar("priorhelpvar","",0.0, 1.0 );
         _addeddummyprior = true;
         RooUniform* priorhelpdist = new RooUniform("priorhelpdist","", *priorhelpvar);
         fPrior = priorhelpdist;
  }

  std::cout << "Imported parameters:\n";
  fParams  = new RooArgList(listPOI);
  const RooArgSet* paramsTmp = params;
  if (paramsTmp!=0)
    fParams->add(*paramsTmp);
  fParams->Print("v");

  fParamsPOI  = new RooArgList(listPOI);
  // create the log-likelihood function
//   fNll = new RooNLLVar("fNll","",*fModel,*fData,true/*extended*/);

  RooArgSet* constrainedParams = fModel->getParameters(*fData);
  fNll = fModel->createNLL(*fData, RooFit::Constrain(*constrainedParams) );

  DefineParameters();

  if(_fillChain){
      SetupRooStatsMarkovChain();
  }
}

// ---------------------------------------------------------
void BCRooInterface::Initialize( const char* rootFile,
             const char* wsName,
             const char* dataName,
             const char* modelName,
             const char* priorName,
             const char* priorNuisanceName,
             const char* paramsName,
             const char* listPOIName )
{
  // retrieve the RooFit inputs from the ROOT file

  /*
  // hard coded names in the workspace
  char* rootFile = "bat_workspace.root";
  char* wsName= "batWS";
  char* dataName= "data";
  char* modelName= "model";
  char* priorName= "priorPOI";
  char* priorNuisanceName= "priorNuisance";
  char* paramsName= "parameters";
  char* listPOIName= "POI";
  */

   std::cout << "Opening " << rootFile << std::endl;
   TFile* file = new TFile(rootFile);
   std::cout << "content :\n";
   file->ls();

   RooWorkspace* bat_ws = (RooWorkspace*) file->Get(wsName);
   bat_ws->Print("v");

   fData = (RooAbsData*) bat_ws->data(dataName);
   fModel = (RooAbsPdf*) bat_ws->function(modelName);

   // make the product of both priors to get the full prior probability function
   RooAbsPdf* priorPOI = (RooAbsPdf*) bat_ws->function(priorName);
   RooAbsPdf* priorNuisance = (RooAbsPdf*) bat_ws->pdf(priorNuisanceName);
   if (priorNuisance!=0 && priorPOI!=0) {
      fPrior = new RooProdPdf("fPrior","complete prior",*priorPOI,*priorNuisance);
   }
   else {
      if ( priorNuisance!=0 )
         fPrior=priorNuisance;
      else if ( priorPOI!=0 )
         fPrior = priorPOI;
      else
         std::cout << "No prior PDF: without taking action the program would crash\n";
         std::cout << "No prior PDF: adding dummy uniform pdf on the interval [0..1]\n";
         priorhelpvar = new RooRealVar("priorhelpvar","",0.0, 1.0 );
         _addeddummyprior = true;
         RooUniform* priorhelpdist = new RooUniform("priorhelpdist","", *priorhelpvar);
         fPrior = priorhelpdist;
   }

   std::cout << "Imported parameters:\n";
   fParams  = new RooArgList(*(bat_ws->set(listPOIName)));
   RooArgSet* paramsTmp = (RooArgSet*) bat_ws->set(paramsName);
   if (paramsTmp!=0) {
      fParams->add(*paramsTmp);
   }
   if (_addeddummyprior == true ) {
      fParams->add(*priorhelpvar);
   }
   fParams->Print("v");

   // create the log-likelihood function
   //fNll = new RooNLLVar("fNll","",*fModel,*fData,true/*extended*/);
   RooArgSet* constrainedParams = fModel->getParameters(*fData);
   fNll = fModel->createNLL(*fData, RooFit::Constrain(*constrainedParams) );

   file->Close();

   DefineParameters();
}


// ---------------------------------------------------------
BCRooInterface::BCRooInterface() : BCModel()
{   // default constructor
   _default_nbins = 500;
   _fillChain = false;
   _addeddummyprior = false;
}

// ---------------------------------------------------------
BCRooInterface::BCRooInterface(const char* name, bool fillChain) : BCModel(name)
{   // another constructor
   _default_nbins = 500;
   _fillChain = fillChain;
}

// ---------------------------------------------------------
BCRooInterface::~BCRooInterface()
{   // default destructor
   if(_fillChain){
      delete _roostatsMarkovChain;
   }
  //test begin
  //cout << "fIterationInterfacecount: " << fIterationInterfacecount << endl;
  //test end
}

// ---------------------------------------------------------
void BCRooInterface::DefineParameters()
{   // define for BAT the list of parameters, range and plot binning

  int nParams = fParams->getSize();
  for (int iParam=0; iParam<nParams; iParam++) {
    RooRealVar* ipar = (RooRealVar*) fParams->at(iParam);
    this->AddParameter(ipar->GetName(),ipar->getMin(),ipar->getMax());
    this->SetNbins(ipar->GetName(),_default_nbins);
    std::cout << "added parameter: " << ipar->GetName() << " defined in range [ " << ipar->getMin() << " - " << ipar->getMax() << " ] with number of bins: " << _default_nbins << " \n";
  }

   for(list< pair<const char*,int> >::iterator listiter = _nbins_list.begin(); listiter != _nbins_list.end(); listiter++){
      this->SetNbins((*listiter).first,(*listiter).second);
      std::cout << "adjusted parameter: " << (*listiter).first << " to number of bins: " << (*listiter).second << " \n";
   }

}

// ---------------------------------------------------------
double BCRooInterface::LogLikelihood(const std::vector<double> & parameters)
{   // this methods returns the logarithm of the conditional probability: p(data|parameters)

   // retrieve the values of the parameters to be tested
  int nParams = fParams->getSize();
  for (int iParam=0; iParam<nParams; iParam++) {
    RooRealVar* ipar = (RooRealVar*) fParams->at(iParam);
    ipar->setVal(parameters.at(iParam));
  }

  // compute the log of the likelihood function
  double logprob = -fNll->getVal();
  return logprob;
}

// ---------------------------------------------------------
double BCRooInterface::LogAPrioriProbability(const std::vector<double> & parameters)
{   // this method returs the logarithm of the prior probability for the parameters: p(parameters).
   // retrieve the values of the parameters to be tested
  int nParams = fParams->getSize();
  for (int iParam=0; iParam<nParams; iParam++) {
    RooRealVar* ipar = (RooRealVar*) fParams->at(iParam);
    ipar->setVal(parameters.at(iParam));
  }
  // compute the log of the prior function
  RooArgSet* tmpArgSet = new RooArgSet(*fParams);
  double prob = fPrior->getVal(tmpArgSet);
  delete tmpArgSet;
  if (prob<1e-300)
    prob = 1e-300;
  return log(prob);
}

// ---------------------------------------------------------
void BCRooInterface::SetNumBins(const char * parname, int nbins)
{
   for(list< pair<const char*,int> >::iterator listiter = _nbins_list.begin(); listiter != _nbins_list.end(); listiter++){
      if(!strcmp((*listiter).first, parname)){
         (*listiter).second = nbins;
         return;
      }
   }
   _nbins_list.push_back( make_pair(parname,nbins) );
}

// ---------------------------------------------------------
void BCRooInterface::SetNumBins(int nbins)
{
   _default_nbins = nbins;
}

// ---------------------------------------------------------
void BCRooInterface::SetupRooStatsMarkovChain()
{
   //RooArgSet * tempRooArgSetPointer = new RooArgSet(* fParams, "paramsMarkovChainPointer");
   //_parametersForMarkovChain = RooArgSet(* fParams, "paramsMarkovChainPointer");
   //fParams->snapshot(_parametersForMarkovChain);
   _parametersForMarkovChainPrevious.add(*fParamsPOI);
   _parametersForMarkovChainCurrent.add(*fParamsPOI);

   cout << "size of _parametersForMarkovChain: " << _parametersForMarkovChainCurrent.getSize() << endl;
   cout << "size of fParamsPOI: " << fParamsPOI->getSize() << endl;
   //cout << "size of tempRooArgSetPointer: " << tempRooArgSetPointer->getSize() << endl;

    _roostatsMarkovChain = new RooStats::MarkovChain();
   //test stuff begin
   //the following way of creating the MArkovChain object does not work!:
   //_roostatsMarkovChain = new RooStats::MarkovChain(_parametersForMarkovChain);
   //test stuff end
   cout << "setting up parameters for RooStats markov chain" << endl;
   _parametersForMarkovChainPrevious.writeToStream(cout, false);

   //initialize fPreviousStep, fCurrentStep, fVecWeights
   int nchains = MCMCGetNChains();
   for(int countChains = 1; countChains<=nchains ; countChains++ ) {
      double tempweight = 1.0;
      fVecWeights.push_back(tempweight);
      vector<double> tempvec;
      TIterator* setiter = fParamsPOI->createIterator();
      double tempval = 0.;
      while(setiter->Next()) {
         tempvec.push_back(tempval);
      }
      fPreviousStep.push_back(tempvec);
      fCurrentStep.push_back(tempvec);
   }

   fFirstComparison = true;

   //test stuff begin
   //var1 = new RooRealVar("var1","var1", 10., 0., 100.);
   //fParamsTest  = new RooArgList();
   //fParamsTest->add(*var1);
   //_parametersForMarkovChain_test.add(*fParamsTest);
   //fIterationInterfacecount = 0;
   //test stuff end
}

// ---------------------------------------------------------
// to be added: add elements with higher weights if elements repeat to avoid unneccessarily long chains
void BCRooInterface::MCMCIterationInterface()
{
  //fIterationInterfacecount+=1;

  if(_fillChain){
      //cout << "Hei-ho running with fillChain!" << endl;
      // get number of chains
      int nchains = MCMCGetNChains();

      // get number of parameters
      int npar = GetNParameters();
      //cout << "number of parameters is: " << npar << endl;

      // loop over all chains and fill histogram
      for (int i = 0; i < nchains; ++i) {
         // get the current values of the parameters. These are
         // stored in fMCMCx.

         TIterator* setiter = fParamsPOI->createIterator();
         int j = 0;

         while(setiter->Next()){

            //check parameter names
            BCParameter * tempBCparam = GetParameter(j);
            //cout << "Chain " << i << "param: " << tempBCparam->GetName() << endl;

            const char * paramnamepointer = (tempBCparam->GetName()).c_str();

            double xij = fMCMCx.at(i * npar + j);
            AddToCurrentChainElement(xij, i, j);
            RooRealVar* parampointer = (RooRealVar*) &(_parametersForMarkovChainCurrent[paramnamepointer]);
            parampointer->setVal(xij);
            j++;
         }
         // will only work if _parametersForMarkovChain had correct info!

//test stuff begin
//var1->setVal( RooRandom::randomGenerator()->Uniform(100.) );
//
//_roostatsMarkovChain->Add(_parametersForMarkovChain_test, 0.001, 1.0);
//
//test stuff end

         if( !(EqualsLastChainElement(i)) ) {
            double weight = GetWeightForChain(i);
            _roostatsMarkovChain->Add(_parametersForMarkovChainPrevious, -1.* MCMCGetLogProbx(j), weight); 
            _parametersForMarkovChainPrevious = _parametersForMarkovChainCurrent;
         }
      }
   }
}

// ---------------------------------------------------------
void BCRooInterface::AddToCurrentChainElement(double xij, int chainNum, int poiNum)
{
   fCurrentStep[chainNum][poiNum] = xij;
}

// ---------------------------------------------------------
bool BCRooInterface::EqualsLastChainElement(int chainNum)
{
   bool equals = true;
   vector<double>::iterator itPrevious = fPreviousStep[chainNum].begin();

   if(fFirstComparison == true) {
      fFirstComparison = false;
      _parametersForMarkovChainPrevious = _parametersForMarkovChainCurrent;
      return true;
   }


   for (vector<double>::iterator itCurrent = fCurrentStep[chainNum].begin(); itCurrent != fCurrentStep[chainNum].end(); ++itCurrent) {
      if(*itCurrent != *itPrevious) {
         equals = false;
         fPreviousStep[chainNum] = fCurrentStep[chainNum];
         break;
      }
      ++itPrevious;
   }

   if(equals == true){
      fVecWeights[chainNum] += 1.0;
   }

   return equals;

}

// ---------------------------------------------------------
double BCRooInterface::GetWeightForChain(int chainNum)
{
   double retval = fVecWeights[chainNum];
   fVecWeights[chainNum]= 1.0 ;
   return retval;
}

