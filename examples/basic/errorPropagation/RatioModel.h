#ifndef __RATIOMODEL__H
#define __RATIOMODEL__H

#include <BAT/BCModel.h>

#include <string>

// ---------------------------------------------------------
class RatioModel : public BCModel
{
public:

    // Constructor and destructor
    RatioModel(const std::string& name);
    ~RatioModel();

    double LogLikelihood(const std::vector<double>& parameters);
    void CalculateObservables(const std::vector<double>& parameters);

};
// ---------------------------------------------------------

#endif

