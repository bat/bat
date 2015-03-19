#ifndef __COMBINATIONMODEL__H
#define __COMBINATIONMODEL__H

#include <BAT/BCModel.h>

// ---------------------------------------------------------
class CombinationModel : public BCModel
{
public:

    // Constructors and destructor
    CombinationModel(const char* name);
    ~CombinationModel();

    // Methods to overload, see file CombinationModel.cxx
    double LogLikelihood(const std::vector<double>& parameters);
};
// ---------------------------------------------------------

#endif

