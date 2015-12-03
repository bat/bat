#include "RatioModel.h"

#include <BAT/BCMath.h>
#include <BAT/BCParameter.h>

#include <cmath>

// ---------------------------------------------------------
RatioModel::RatioModel(const std::string& name)
    : BCModel(name)
{
    // define the parameters x and y
    AddParameter("x", 0., 8.); // index 0
    AddParameter("y", 0., 16.); // index 1

    GetParameters().SetPriorConstantAll();

    AddObservable("r", 0, 2, "#frac{x}{y}");
}

// ---------------------------------------------------------
RatioModel::~RatioModel()
{
}

// ---------------------------------------------------------
double RatioModel::LogLikelihood(const std::vector<double>& parameters)
{
    // model = Gaus(x | mean = 4, std. dev. = 1) * Gaus(y | mean = 8, std. dev. = 2)
    return BCMath::LogGaus(parameters[0], 4, 1) + BCMath::LogGaus(parameters[1], 8, 2);
}

// ---------------------------------------------------------
void RatioModel::CalculateObservables(const std::vector<double>& parameters)
{
    // store ratio, if demoninator is not zero
    if (parameters[1] != 0)
        GetObservable(0).Value(parameters[0] / parameters[1]);
    // else store zero
    else if (parameters[0] == 0)
        GetObservable(0).Value(0);
}

