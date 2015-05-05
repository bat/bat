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

    PolAsymm(std::string name);

    ~PolAsymm()
    { /* empty destructor */ }

    // necessary to overload pure virtual BCFitter function
    bool Fit()
    { return false; }

    // necessary to overload pure virtual BCFitter function
    void DrawFit(const char* options, bool flaglegend = false)
    { }

    // fit function returning expectation value for each data point
    double FitFunction(const std::vector<double>& x, const std::vector<double>& par);

    // loglikelihood function - probability of the data given the parameters
    double LogLikelihood(const std::vector<double>& par);

};

#endif

