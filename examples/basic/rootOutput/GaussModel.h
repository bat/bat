#ifndef __GAUSSMODEL__H
#define __GAUSSMODEL__H

#include <BAT/BCModel.h>

#include <string>

class GaussModel : public BCModel
{
public:

    GaussModel(const std::string& name);

    ~GaussModel()
    { /* empty destructor */ }

    // Method to overload, see file GaussModel.cxx
    double LogLikelihood(const std::vector<double>& parameters);
};
// ---------------------------------------------------------

#endif

