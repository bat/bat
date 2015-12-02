#ifndef __BAT__REFERENCECOUNTING__H
#define __BAT__REFERENCECOUNTING__H

#include <BAT/BCModel.h>

#include <string>
#include <vector>

class ReferenceCounting : public BCModel
{
public:

    // An enumerator for the prior evaluation options.
    enum EvalOption {
        kAnalytic,              // use analytic approximation
        kHistogram,             // fill a histogram with the prior before running
        kApprox                 // use approximation
    };

    // Constructor
    ReferenceCounting(const std::string& name, unsigned nobs = 0);

    // Destructor
    ~ReferenceCounting()
    { /* empty destructor */ }

    double LogLikelihood(const std::vector<double>& parameters);

    // set option of how to evaluate the reference prior
    bool SetPrior(ReferenceCounting::EvalOption option, double shape, double rate);

    void SetNObs(unsigned n)
    { fNObs = n; }

    unsigned GetNObs()
    { return fNObs; }

protected:

    unsigned fNObs;                  // number of observed events

};
// ---------------------------------------------------------

#endif

