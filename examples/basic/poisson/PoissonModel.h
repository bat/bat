#ifndef __POISSONMODEL__H
#define __POISSONMODEL__H

#include <BAT/BCModel.h>

#include <string>

class PoissonModel : public BCModel
{
public:

    PoissonModel(const std::string& name);

    ~PoissonModel()
    { /* empty desctructor */ }

    // set number of observed events
    void SetNObs(unsigned nobs);

    // get number of observed events
    unsigned GetNObs() const
    { return fNObs; }

    double LogLikelihood(const std::vector<double>& parameters);

private:

    // number of observed events
    unsigned fNObs;
};
// ---------------------------------------------------------

#endif

