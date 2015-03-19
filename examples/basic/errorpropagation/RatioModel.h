#ifndef __RATIOMODEL__H
#define __RATIOMODEL__H

#include <BAT/BCModel.h>

// ---------------------------------------------------------
class RatioModel : public BCModel
{
public:

    // Constructor and destructor
    RatioModel(const char* name);
    ~RatioModel();

    double LogLikelihood(const std::vector<double>& parameters);
    void CalculateObservables(const std::vector<double>& parameters);

};
// ---------------------------------------------------------

#endif

