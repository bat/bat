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
		     const char* wsName = "bat_ws",
		     const char* dataName = "fData",
		     const char* modelName = "fModel",
		     const char* priorName = "fPrior",
		     const char* observablesName = "fObservables",
		     const char* paramsName = "fParams" );

	private:

		RooDataSet* fData;        // data to test
		RooAbsPdf*  fModel;       // likelihood model describing the observables
		RooNLLVar*  fNll;         // pointer to negative log-likelihood function
		RooArgSet*  fObservables; // list of observables measured for each event
		RooArgList* fParams;      // list of parameters of interest
		RooAbsPdf*  fPrior;       // function describing the prior probability of the parameters

};

#endif

