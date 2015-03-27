
#ifndef __BC_TEST__GAUSSMODEL__H
#define __BC_TEST__GAUSSMODEL__H

#include <BAT/BCModel.h>

// ---------------------------------------------------------
class GaussModel : public BCModel
{
public:

    // Constructors and destructor
    GaussModel(const char* name, const unsigned& nParameters, unsigned long loopIterations = 0);
    virtual ~GaussModel();

    // Methods to overload, see file GaussModel.cxx
    virtual double LogLikelihood(const std::vector<double>& parameters);

    unsigned long Calls() const
    {
        return fCalls;
    }

    void MCMCCurrentPointInterface(std::vector<double>& /*p*/, int /*c*/, bool /*accepted*/);

private:
    /**
     * Used in likelihood to prolong artificially.
     */
    unsigned long fLoopIterations;

    /**
     * Count how often likelihood is called
     */
    unsigned long fCalls;

};
// ---------------------------------------------------------

#endif