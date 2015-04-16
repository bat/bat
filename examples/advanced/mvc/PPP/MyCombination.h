#ifndef __MYCOMBINATION__H
#define __MYCOMBINATION__H

#include <BAT/BCMVCombination.h>

#include <vector>

// ---------------------------------------------------------
class MyCombination : public BCMVCombination
{
public:

    // Constructor
    MyCombination();

    // Destructor
    ~MyCombination();

    void SetFlagPhysicalConstraints(bool flag)
    { fFlagPhysicalConstraints = flag; };

    // BAT methods

    double LogLikelihood(const std::vector<double>& parameters);

    void CalculateObservables(const std::vector<double>& pars);

private:

    // flag for imposing physical constraints or not
    bool fFlagPhysicalConstraints;

};
// ---------------------------------------------------------

#endif

