/*
 * Copyright (C) 2007-2015, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

#include "GaussModel.h"

#include <test.h>
#include <BAT/BCMath.h>
#include <BAT/BCParameter.h>

using namespace test;

// ---------------------------------------------------------
GaussModel::GaussModel(const std::string& name, const unsigned& nParameters, unsigned long loopIterations) :
    BCModel(name),
    fLoopIterations(loopIterations),
    fCalls(0)
{
    // add identical, independent parameters
    for (unsigned i = 0; i < nParameters ; ++i) {
        std::string parName("par");
        AddParameter((parName + test::stringify(i)).c_str(), -15.0, 15.0);
    }
    SetPriorConstantAll();
}

// ---------------------------------------------------------
GaussModel::~GaussModel()
{
}

void GaussModel::MCMCCurrentPointInterface(const std::vector<double>& /*p*/, int /*c*/, bool /*accepted*/)
{
    #pragma omp critical(GaussModel__LogLikelihood)
    {
        ++fCalls;
    }
}

// ---------------------------------------------------------
double GaussModel::LogLikelihood(const std::vector<double>& parameters)
{
    // run extra loop to make likelihood evaluation slower(!)
    if (fLoopIterations) {
        for (unsigned long i = 0 ; i < fLoopIterations ; ++i)
            BCMath::LogFact(10);
    }

    // check that fixed parameters are indeed fixed to the right value
    for (unsigned i = 0; i < parameters.size(); i++)
        if (fParameters[i].Fixed())
            TEST_CHECK_NEARLY_EQUAL(fParameters[i].GetFixedValue(), parameters[i], 1e-15);

    // assume a normalized Gaussian Likelihood with N independent variables
    static const double normalized = true;
    double logprob = 0;
    for (unsigned i = 0; i < parameters.size(); i++)
        logprob += BCMath::LogGaus(parameters[i], mean(), sigma(), normalized);

    return logprob;
}
