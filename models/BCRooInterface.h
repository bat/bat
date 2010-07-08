#ifndef __BCROOINTERFACE__H
#define __BCROOINTERFACE__H

#include <BAT/BCModel.h>

class RooAbsReal;
class RooAbsData;
class RooAbsPdf;
class RooArgSet;
class RooArgList;

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
		void Initialize( RooAbsData& data,
				 RooAbsPdf& model,
				 RooAbsPdf& prior,
				 RooAbsPdf* priorNuisance,
				 RooArgSet* params,
				 RooArgSet& listPOI );

	private:

		RooAbsData* fData;        // data to test
		RooAbsPdf*  fModel;       // likelihood model describing the observables
		RooAbsReal*  fNll;         // pointer to negative log-likelihood function
		RooArgSet*  fObservables; // list of observables measured for each event
		RooArgList* fParams;      // list of parameters of interest
		RooAbsPdf*  fPrior;       // function describing the prior probability of the parameters

};

#endif
