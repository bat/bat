
#ifndef __BC_TEST__GAUSSMODEL__H
#define __BC_TEST__GAUSSMODEL__H

#include <BAT/BCModel.h>

/**
 * @class GaussModel A unit Gaussian.
 */
// ---------------------------------------------------------
class GaussModel : public BCModel
{
public:

    // Constructors and destructor
    GaussModel(const std::string& name, const unsigned& nParameters, unsigned long loopIterations = 0);
    virtual ~GaussModel();

    // Methods to overload, see file GaussModel.cxx
    virtual double LogLikelihood(const std::vector<double>& parameters);

    unsigned long Calls() const
    {
        return fCalls;
    }

    double mean() const
    {
        return 0;
    }

    double sigma() const
    {
        return 1;
    }

    void MCMCCurrentPointInterface(const std::vector<double>& /*p*/, int /*c*/, bool /*accepted*/);

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
