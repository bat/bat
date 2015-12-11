/*!
 * \class BAT::PerfTestMCMC
 * \brief A performance test class for BAT
 */

/*
 * Copyright (C) 2009, Daniel Kollar and Kevin Kroeninger.
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */


#ifndef BAT_PERFTESTMCMC
#define BAT_PERFTESTMCMC

#include <include/PerfTest.h>

#include <BAT/BCModel.h>

#include <TF1.h>
#include <TGraph.h>
#include <TH2.h>

#include <string>
#include <vector>

class PerfTestMCMC : public PerfTest, public BCModel
{

public:

    /** \name Constructors and destructors  */
    /* @{ */

    /** The default constructor */
    PerfTestMCMC(const std::string& name = "unknown");

    /** The default destructor */
    ~PerfTestMCMC();

    /* @} */

    /** Set the variation parameter.
     * @param par the parameter value
     * @param name the name of the varied parameter.
     * @return an error code. */
    virtual int SetVarPar(double value, const std::string& name);

    virtual void SetProposal(bool multivariate, double dof)
    {
        SetProposeMultivariate(multivariate);
        SetProposalFunctionDof(dof);
    }

    /** Run before test.
     * @return an error code. */
    int PreTest();

    /** Run after test.
     * @return an error code. */
    int PostTest();

    /** Run the test.
     * @return an error code. */
    int RunTest();

    /** Defines the subtests. */
    void DefineSubtests();

    /** Writes the test to file.
     * @return an error code. */
    int WriteResults();

    /** Define precision settings. */
    void PrecisionSettings(PerfTest::Precision);

    /* @} */

    // inherited methods
    void MCMCUserIterationInterface();

private:

    std::vector<TGraph*> fCorrelation;
    std::vector<TH2D*> fHistCorr;
    std::vector<std::vector<double> > fXOld;
};

#endif
