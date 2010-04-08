
#include <RooGlobalFunc.h>
#include <RooMsgService.h>
#include <RooProdPdf.h>
#include <RooRealVar.h>
#include <RooWorkspace.h>

#include <RooAbsData.h>
#include <RooAbsPdf.h>
#include <RooNLLVar.h>
#include <RooArgSet.h>
#include <RooArgList.h>

#include <TFile.h>

#include <iostream>

#include <BAT/BCMath.h>

#include "BCRooInterface.h"

// ---------------------------------------------------------
void BCRooInterface::Initialize(
		const char * rootFile,
		const char * wsName,
		const char * dataName,
		const char * modelName,
		const char * priorName,
		const char * priorNuisanceName,
		const char * paramsName,
		const char * listPOIName )
{
  // retrieve the RooFit inputs from the ROOT file

  std::cout << "Opening " << rootFile << std::endl;
  TFile * file = new TFile(rootFile);
  std::cout << "content :\n";
  file->ls();
  
  RooWorkspace * bat_ws = (RooWorkspace *) file->Get(wsName);
  bat_ws->Print("v");
  
  fData = (RooAbsData *) bat_ws->data(dataName);
  fModel = (RooAbsPdf *) bat_ws->function(modelName);
  
  // make the product of both priors to get the full prior probability function
  RooAbsPdf * priorPOI = (RooAbsPdf *) bat_ws->function(priorName);
  RooAbsPdf * priorNuisance = (RooAbsPdf *) bat_ws->pdf(priorNuisanceName);
  if (priorNuisance!=0 && priorPOI!=0) {
    fPrior = new RooProdPdf("fPrior","complete prior",*priorPOI,*priorNuisance);
  }
  else {
    if ( priorNuisance!=0 )
	 	fPrior=priorNuisance;
    else if ( priorPOI!=0 )
	 	fPrior = priorPOI;
    else
	 	std::cout << "No prior PDF: the program will crash\n";
  }
  
  std::cout << "Imported parameters:\n";
  fParams  = new RooArgList(*(bat_ws->set(listPOIName)));
  RooArgSet * paramsTmp = (RooArgSet *) bat_ws->set(paramsName);
  if (paramsTmp!=0)
    fParams->add(*paramsTmp);
  fParams->Print("v");

  // create the log-likelihood function
  fNll = new RooNLLVar("fNll","",*fModel,*fData,true/*extended*/);
  
  file->Close();
  
  DefineParameters();
}

// ---------------------------------------------------------
BCRooInterface::BCRooInterface() : BCModel()
{	// default constructor

}

// ---------------------------------------------------------
BCRooInterface::BCRooInterface(const char * name) : BCModel(name)
{	// another constructor

}

// ---------------------------------------------------------
BCRooInterface::~BCRooInterface()
{	// default destructor

}

// ---------------------------------------------------------
void BCRooInterface::DefineParameters()
{	// define for BAT the list of parameters, range and plot binning

  int default_nbins = 500;
	
  int nParams = fParams->getSize();
  for (int iParam=0; iParam<nParams; iParam++) {
    RooRealVar * ipar = (RooRealVar *) fParams->at(iParam);
    this->AddParameter(ipar->GetName(),ipar->getMin(),ipar->getMax());
    this->SetNbins(ipar->GetName(),default_nbins);
    std::cout << "added parameter: " << ipar->GetName() << " defined in range [ " << ipar->getMin() << " - " << ipar->getMax() << " ]\n";
  }
}

// ---------------------------------------------------------
double BCRooInterface::LogLikelihood(std::vector <double> parameters)
{	// this methods returns the logarithm of the conditional probability: p(data|parameters)

	// retrieve the values of the parameters to be tested
  int nParams = fParams->getSize();
  for (int iParam=0; iParam<nParams; iParam++) {
    RooRealVar * ipar = (RooRealVar *) fParams->at(iParam);
    ipar->setVal(parameters.at(iParam));
  }

  // compute the log of the likelihood function
  double logprob = -fNll->getVal();
  return logprob;
}

// ---------------------------------------------------------
double BCRooInterface::LogAPrioriProbability(std::vector <double> parameters)
{	// this method returs the logarithm of the prior probability for the parameters: p(parameters).

	// retrieve the values of the parameters to be tested
  int nParams = fParams->getSize();
  for (int iParam=0; iParam<nParams; iParam++) {
    RooRealVar * ipar = (RooRealVar *) fParams->at(iParam);
    ipar->setVal(parameters.at(iParam));
  }

  // compute the log of the prior function
  RooArgSet * tmpArgSet = new RooArgSet(*fParams);
  double prob = fPrior->getVal(tmpArgSet);
  delete tmpArgSet;
  if (prob<1e-300)
    prob = 1e-300;
  return log(prob);
}

