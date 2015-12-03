#ifndef __COMBINATIONMODEL__H
#define __COMBINATIONMODEL__H

#include <BAT/BCModel.h>

#include <string>

class CombinationModel : public BCModel
{
public:

    CombinationModel(const std::string& name,
                     double new_mean, double new_sigma,
                     double old_mean, double old_sigma);

    ~CombinationModel()
    { /* empty destructor */ }

    double LogLikelihood(const std::vector<double>& pars);

protected:

    double fNewMean;            // newly measured mean
    double fNewSigma;           // newly measured sigma
};

#endif

