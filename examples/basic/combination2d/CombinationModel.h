#ifndef __COMBINATIONMODEL__H
#define __COMBINATIONMODEL__H

#include <BAT/BCModel.h>

#include <string>

class CombinationModel : public BCModel
{
public:

    CombinationModel(const std::string& name,
                     double new_mass_mean, double new_mass_sigma,
                     double new_xs_mean, double new_xs_sigma,
                     double rho,
                     double old_mass_mean, double old_mass_sigma,
                     double old_xs_mean, double old_xs_sigma);

    ~CombinationModel()
    {/* empty destructor*/ }

    double LogLikelihood(const std::vector<double>& parameters);

protected:

    double fMassMean;           /* measured mass mean */
    double fMassSigma;          /* measured mass std. dev. */
    double fXSMean;             /* measured xs mean */
    double fXSSigma;            /* measured xs std. dev. */
    double fRho;                /* correlation */

};

#endif

