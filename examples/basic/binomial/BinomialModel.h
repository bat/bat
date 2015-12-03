#ifndef __BINOMIALMODEL__H
#define __BINOMIALMODEL__H

#include <BAT/BCModel.h>

// ---------------------------------------------------------
class BinomialModel : public BCModel
{
public:

    BinomialModel(const std::string& name, unsigned ntotal, unsigned nselected);

    ~BinomialModel();

    double LogLikelihood(const std::vector<double>& pars);

protected:

    unsigned fNTotal;          // total number of events
    unsigned fNSelected;       // selected (observed) number of events
};
// ---------------------------------------------------------

#endif

