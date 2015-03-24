/*
 * Copyright (C) 2007-2015, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

#include "test.h"
#include "GaussModel.h"

using namespace test;

namespace
{
GaussModel factory(const char* name = "BCEngineMCMC_TEST-model")
{
    GaussModel m(name, 5);

    m.MCMCSetNChains(2);
    m.MCMCSetNIterationsEfficiencyCheck(500);
    m.MCMCSetNIterationsConvergenceCheck(500);
    m.MCMCSetNIterationsClearConvergenceStats(20000);
    m.MCMCSetNIterationsPreRunMax(1000000);
    m.MCMCSetNIterationsPreRunMin(1000);
    m.MCMCSetNIterationsRun(10000);
    m.MCMCSetMinimumEfficiency(0.15);
    m.MCMCSetMaximumEfficiency(0.35);

    return m;
}
}

class MultivariateProposalTest :
    public TestCase
{
public:
    MultivariateProposalTest() :
        TestCase("MultivariateProposal test")
    {
        BCLog::SetLogLevelScreen(BCLog::debug);
    }

    virtual void run() const
    {
        TEST_SECTION(
            "convergence checking",

            GaussModel m = ::factory();
            m.MCMCSetCorrectRValueForSamplingVariability(true);
            m.MCMCSetMultivariateProposalFunction(true);

            m.MCMCSetRandomSeed(915);

            m.MarginalizeAll(BCIntegrate::kMargMetropolis);

            // no accidental calls to likelihood
            TEST_CHECK_EQUAL(m.Calls(), m.MCMCGetNChains() * (m.MCMCGetNIterationsConvergenceGlobal() + m.MCMCGetNIterationsRun()));
        );
    }
} multivariateProposalTest;

// Local Variables:
// compile-command: "make -C .. check TESTS= && (./BCEngineMCMC_TEST || cat test-suite.log)"
// End:
