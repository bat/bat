#ifndef __BCROOINTERFACE__H
#define __BCROOINTERFACE__H

#include "RooAbsPdf.h"
#include "RooArgList.h"
#include "RooArgSet.h"
#include "RooDataSet.h"
#include "RooNLLVar.h"

#include <BAT/BCModel.h>

/* ---------------------------------------------------------
 * interface allowing to run BAT on a problem/data defined in a
 * standard RooFit workspace format
 * --------------------------------------------------------- */

class BCRooInterface : public BCModel
{

	public:

		// Constructors and destructor
		BCRooInterface( );

		BCRooInterface( const char* name );

		~BCRooInterface();

		// Overloaded methods
		void DefineParameters();
		double LogAPrioriProbability(std::vector <double> parameters);
		double LogLikelihood(std::vector <double> parameters);

		// Other method of this class
		void Initialize( const char* rootFile,
		     const char* wsName = "batWS",
		     const char* dataName = "data",
		     const char* modelName = "model",
		     const char* priorName = "priorPOI",
		     const char* priorNuisanceName= "priorNuisance",
		     const char* paramsName = "parameters",
		     const char* listPOIName = "POI" );

	private:

		RooAbsData* fData;        // data to test
		RooAbsPdf*  fModel;       // likelihood model describing the observables
		RooNLLVar*  fNll;         // pointer to negative log-likelihood function
		RooArgSet*  fObservables; // list of observables measured for each event
		RooArgList* fParams;      // list of parameters of interest
		RooAbsPdf*  fPrior;       // function describing the prior probability of the parameters

};

#endif

