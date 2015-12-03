#include "CombinationModel.h"

#include <BAT/BCGaussianPrior.h>
#include <BAT/BCMath.h>
#include <BAT/BCPositiveDefinitePrior.h>

// ---------------------------------------------------------
CombinationModel::CombinationModel(const std::string& name,
                                   double new_mass_mean, double new_mass_sigma,
                                   double new_xs_mean, double new_xs_sigma,
                                   double rho,
                                   double old_mass_mean, double old_mass_sigma,
                                   double old_xs_mean, double old_xs_sigma)
    : BCModel(name),
      fMassMean(new_mass_mean),
      fMassSigma(new_mass_sigma),
      fXSMean(new_xs_mean),
      fXSSigma(new_xs_sigma),
      fRho(rho)
{
    // add parameter for mass, with range large enough to reach out to
    // 4 sigma in either direction from the old and new means
    // LaTeX name "Mass" and units "MeV"
    AddParameter("mass",
                 std::min<double>(fMassMean - 4 * fMassSigma, old_mass_mean - 4 * old_mass_sigma),
                 std::max<double>(fMassMean + 4 * fMassSigma, old_mass_mean + 4 * old_mass_sigma),
                 "Mass", "[MeV]");

    // add parameter for cross section, with range large enough to
    // reach out to 4 sigma in either direction from the old and new means
    // LaTeX name "#sigma" and units "ab"
    AddParameter("xs",
                 std::min<double>(fXSMean - 4 * fXSSigma, old_xs_mean - 4 * old_xs_sigma),
                 std::max<double>(fXSMean + 4 * fXSSigma, old_xs_mean + 4 * old_xs_sigma),
                 "#sigma", "[ab]");

    // set priors
    GetParameter("mass").SetPrior(new BCPositiveDefinitePrior(new BCGaussianPrior(old_mass_mean, old_mass_sigma)));
    GetParameter("xs").SetPrior(new BCPositiveDefinitePrior(new BCGaussianPrior(old_xs_mean, old_xs_sigma)));
};

// ---------------------------------------------------------
double CombinationModel::LogLikelihood(const std::vector<double>& pars)
{
    // bivariate gaussian with correlation fRho
    return BCMath::LogGaus(pars[0], fMassMean, fMassMean * sqrt(1 - fRho * fRho))
           + BCMath::LogGaus(pars[1], fXSMean, fXSSigma * sqrt(1 - fRho * fRho))
           + fRho / (1 - fRho * fRho) * (pars[0] - fMassMean) / fMassMean * (pars[1] - fXSMean) / fXSSigma
           - 0.5 * log(1 - fRho * fRho);
}
