/*
 * Copyright (C) 2007-2018, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

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
    virtual void CalculateObservables(const std::vector<double>& parameters);

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

    double evidence() const;

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
