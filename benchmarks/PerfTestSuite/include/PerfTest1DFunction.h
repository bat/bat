/*!
 * \class BAT::PerfTest1DFunction
 * \brief A performance test class for BAT
 */

/*
 * Copyright (C) 2009, Daniel Kollar and Kevin Kroeninger.
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */


#ifndef BAT_PERFTEST1DFUNCTION
#define BAT_PERFTEST1DFUNCTION

#include <include/PerfTestMCMC.h>
#include <include/PerfTest.h>

#include <TF1.h>
#include <TH1.h>

#include <math.h>
#include <string>
#include <vector>


class PerfTest1DFunction : public PerfTestMCMC
{

public:

    /** \name Constructors and destructors  */
    /* @{ */

    /** The default constructor */
    PerfTest1DFunction(const std::string& name = "unknown", TF1* func = 0);

    /** The default destructor */
    virtual ~PerfTest1DFunction();

    /* @} */

    /** Run after the test
     * @return an error code. */
    int PostTest();

    /** Defines the subtests. */
    void DefineSubtests();

    /* @} */

    // inherited methods
    double LogAPrioriProbability(const std::vector<double>& /*pars*/)
    { return 0; }

    double LogLikelihood(const std::vector<double>& pars)
    { return log(fFunction->Eval(pars[0])); }

private:

    /** The test function. */
    TF1* fFunction;
};

#endif
