#ifndef __POL2ASYMM__H
#define __POL2ASYMM__H

/*
 * This class derives from Pol1Asymm. It describes a linear
 * correlation relation between measured points. Two parameters
 * are defined within the model, an offset and a slope.
 * The data are points (x,y) with an asymmetric uncertainty on y.
 * The uncertainty is assumed to be two half gaussians with different
 * widths.
 */

// ---------------------------------------------------------

#include <vector>

#include "Pol1Asymm.h"

// ---------------------------------------------------------

class Pol2Asymm : public Pol1Asymm
{
public:
        // default constructor
        Pol2Asymm();
        // constructor setting the name of the model
        Pol2Asymm(const char * name);
        // destructor
        ~Pol2Asymm();

        // define parameters of the model
        virtual void DefineParameters();

        // fit function returning expectation value for each data point
        virtual double FitFunction(const std::vector<double> &x, const std::vector<double> & par);

        // NOTE ON LogLikelihood:
        // ----------------------
        // We want this model to have the same definition of uncertainties
        // as the Pol1Asymm so no new LogLikelihood definition is needed.
        // The one in Pol1Asymm is general enough. The only thing new here
        // is the FitFunction which we have defined.

        // NOTE ON LogAPrioriProbability:
        // ----------------------
        // If we want to use flat prior in all the variables like in
        // Pol1Asymm we don't have to define LogAPrioriProbability again
        // since it's definition in Pol1Asymm model was general enough
};

// ---------------------------------------------------------

#endif

