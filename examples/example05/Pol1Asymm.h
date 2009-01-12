#ifndef __POL1ASYMM__H
#define __POL1ASYMM__H

/*
 * This class derives from BCModel. It describes a linear
 * correlation relation between measured points. Two parameters
 * are defined within the model, an offset and a slope.
 * The data are points (x,y) with an asymmetric uncertainty on y.
 * The uncertainty is assumed to be two half gaussians with different
 * widths.
 */

// ---------------------------------------------------------

#include <vector>

#include <BAT/BCModel.h>

// ---------------------------------------------------------

class Pol1Asymm : public BCModel
{
	public:
		// default constructor
		Pol1Asymm();
		// constructor setting the name of the model
		Pol1Asymm(const char * name);
		// destructor
		~Pol1Asymm();

		// define parameters of the model
		virtual void DefineParameters();

		// fit function returning expectation value for each data point
		virtual double FitFunction(std::vector <double> x, std::vector <double> par);

		// loglikelihood function - probability of the data given parameters
		virtual double LogLikelihood(std::vector <double> par);

		// prior probability
		virtual double LogAPrioriProbability(std::vector <double> par);
};

// All methods were defined virtual so that we can derive from this class
// and overload them as needed if necessary. This has been done in
// class Pol2Asymm which is derived from this one.

// ---------------------------------------------------------

#endif

