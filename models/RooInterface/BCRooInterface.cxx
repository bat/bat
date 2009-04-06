#include "RooGlobalFunc.h"
#include "RooMsgService.h"
#include "RooRealVar.h"
#include "RooWorkspace.h"
#include "TFile.h"

#include "BCRooInterface.h"
#include <BAT/BCMath.h>
#include <iostream>


// ---------------------------------------------------------
void BCRooInterface::Initialize( const char* rootFile,
		const char* wsName,
		const char* dataName,
		const char* modelName,
		const char* priorName,
		const char* observablesName,
		const char* paramsName )
{	// retrieve the RooFit inputs from the ROOT file

	TFile* file = new TFile(rootFile);
	RooWorkspace* bat_ws = (RooWorkspace*) file->Get(wsName);
	bat_ws->Print("v");

	fData = (RooDataSet*) bat_ws->data(dataName);
	fModel = (RooAbsPdf*) bat_ws->function(modelName);
	fPrior = (RooAbsPdf*) bat_ws->function(priorName);

	// START: temporary fix until RooWorkspace supports RooArgList and RooArgSet import
	RooArgSet* observablesTmp = (RooArgSet*) file->Get(observablesName);
	RooArgList observablesTmp2(*observablesTmp);
	fObservables = new RooArgSet();
	int nParams = observablesTmp2.getSize();
	for (int iParam=0; iParam<nParams; iParam++)
	{
		fObservables->add(*((RooRealVar*)bat_ws->fundArg(((RooRealVar*) observablesTmp2.at(iParam))->GetName())));
	}

	RooArgList* paramsTmp = (RooArgList*) file->Get(paramsName);
	fParams = new RooArgList();
	nParams = paramsTmp->getSize();
	for (int iParam=0; iParam<nParams; iParam++)
	{
		fParams->add(*((RooRealVar*)bat_ws->fundArg(((RooRealVar*) paramsTmp->at(iParam))->GetName())));
	}
	// END: temporary fix until RooWorkspace supports RooArgList and RooArgSet import

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
BCRooInterface::BCRooInterface(const char* name) : BCModel(name)
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
		RooRealVar* ipar = (RooRealVar*) fParams->at(iParam);
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
		RooRealVar* ipar = (RooRealVar*) fParams->at(iParam);
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
		RooRealVar* ipar = (RooRealVar*) fParams->at(iParam);
		ipar->setVal(parameters.at(iParam));
	}

	// compute the log of the prior function
	RooArgSet* tmpArgSet = new RooArgSet(*fParams);
	double prob = fPrior->getVal(tmpArgSet);
	delete tmpArgSet;
	if (prob<1e-300) prob = 1e-300;
	return log(prob);
}

