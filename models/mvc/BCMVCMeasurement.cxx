#include "BCMVCMeasurement.h"

#include <cmath>

// ---------------------------------------------------------
BCMVCMeasurement::BCMVCMeasurement(std::string name) : fName(name)
						     , fObservable(-1)
						     , fCentralValue(0.)
						     , fUncertainties(0)
						     , fFlagActive(true)
{
}

// ---------------------------------------------------------
BCMVCMeasurement::~BCMVCMeasurement()
{
}

// ---------------------------------------------------------
double BCMVCMeasurement::GetTotalUncertainty()
{
  int n = int(fUncertainties.size());

  double sum2 = 0;

  // sum all uncertainties in quadrature
  for (int i = 0; i < n; ++i) {
    double unc = GetUncertainty(i);
    sum2+=unc*unc;
  }

  return sqrt(sum2);
}

// ---------------------------------------------------------
