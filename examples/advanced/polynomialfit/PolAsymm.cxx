#include "PolAsymm.h"

#include <BAT/BCDataPoint.h>
#include <BAT/BCDataSet.h>
#include <BAT/BCLog.h>
#include <BAT/BCMath.h>

#include <TString.h>

#include <cmath>

// ---------------------------------------------------------
PolAsymm::PolAsymm(std::string name)
    : BCFitter(name)
{
}

// ---------------------------------------------------------
double PolAsymm::FitFunction(const std::vector<double>& x, const std::vector<double>& par)
{
    // n-th order polynomial
    double r = 0;
    double X = 1;
    for (unsigned i = 0; i < par.size(); ++i) {
        r += par[i] * X;
        X *= x[0];
    }

    return r;
}

// ---------------------------------------------------------
double PolAsymm::LogLikelihood(const std::vector<double>& par)
{
    double logl = 0.;

    // loop over the data points
    for (unsigned i = 0 ; i < GetNDataPoints(); i++) {
        // get data point
        std::vector<double> x = GetDataSet()->GetDataPoint(i)->GetValues();

        // calculate the value of the function at this point
        double y_func = FitFunction(x, par);

        // Likelihood *= asymmetric Gaussian evaluated at y_func given
        // mode from the data point (x[1]) and standard deviation
        // above mode (x[2]) and below mode (x[3]); including normalization constant
        logl += BCMath::LogSplitGaus(y_func, x[1], x[2], x[3], true);
    }

    return logl;
}
