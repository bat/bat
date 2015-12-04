#ifndef __POLASYMM__H
#define __POLASYMM__H

/*
 * This class derives from BCModel. It describes a polynomial
 * relationship between measured points. The data are points (x,y)
 * with an asymmetric uncertainty on y (above and below y).
 */

#include "BAT/BCFitter.h"

#include <string>
#include <vector>

class PolAsymm : public BCFitter
{
public:

    /** Create a polynomial model with name and the parameters. */
    PolAsymm(const std::string& name, const BCParameterSet& parameters);

    virtual ~PolAsymm()
    { /* empty destructor */ }

    // necessary to overload pure virtual BCFitter function
    virtual void Fit() {}

    // necessary to overload pure virtual BCFitter function
    virtual void DrawFit(const std::string& options, bool flaglegend = false)
    {}

    // fit function returning expectation value for each data point
    virtual double FitFunction(double* x, double* par);

    // loglikelihood function - probability of the data given the parameters
    double LogLikelihood(const std::vector<double>& par);
};

#endif
