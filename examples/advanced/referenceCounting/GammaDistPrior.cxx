#include "GammaDistPrior.h"

#include <TMath.h>

// ---------------------------------------------------------
GammaDistPrior::GammaDistPrior(double shape, double rate)
    : BCPrior(),
      fShape(shape),
      fRate(rate)
{
}

// ---------------------------------------------------------
GammaDistPrior::GammaDistPrior(const GammaDistPrior& other)
    : BCPrior(other),
      fShape(other.fShape),
      fRate(other.fRate)
{
}

// ---------------------------------------------------------
double GammaDistPrior::GetLogPrior(double x)
{
    return fShape * log(fRate) + (fShape - 1) * log(x) - fRate * x - TMath::LnGamma(fShape);
}
