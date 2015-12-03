#include "CombinationModel.h"

#include <BAT/BCGaussianPrior.h>
#include <BAT/BCMath.h>
#include <BAT/BCPositiveDefinitePrior.h>

// ---------------------------------------------------------
CombinationModel::CombinationModel(const std::string& name,
                                   double new_mean, double new_sigma,
                                   double old_mean, double old_sigma)
    : BCModel(name),
      fNewMean(new_mean),
      fNewSigma(new_sigma)
{
    // add parameter for mass, with range large enough to reach out to
    // 4 sigma in either direction from the old and new means
    // LaTeX name "mass" and units "MeV"
    AddParameter("mass",
                 std::min<double>(fNewMean - 4 * fNewSigma, old_mean - 4 * old_sigma),
                 std::max<double>(fNewMean + 4 * fNewSigma, old_mean + 4 * old_sigma),
                 "mass", "[MeV]");

    // set its prior to the old mean and sigma
    GetParameter("mass").SetPrior(new BCPositiveDefinitePrior(new BCGaussianPrior(old_mean, old_sigma)));
}

// ---------------------------------------------------------
double CombinationModel::LogLikelihood(const std::vector<double>& pars)
{
    return BCMath::LogGaus(pars[0], fNewMean, fNewSigma);
}
