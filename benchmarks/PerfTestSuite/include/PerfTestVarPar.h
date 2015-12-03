/*!
 * \class BAT::PerfTestVarPar
 * \brief A performance test class for BAT
 */

/*
 * Copyright (C) 2009, Daniel Kollar and Kevin Kroeninger.
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */


#ifndef BAT_PERFTESTVARPAR
#define BAT_PERFTESTVARPAR

#include <string>
#include <vector>

#include <include/PerfTest.h>
#include <include/PerfTestMCMC.h>

class TH1D;
class TCanvas;
class TGraphErrors;

class PerfTestVarPar : public PerfTest
{

public:

    /** \name Enumerators  */
    /* @{ */

    /* @} */
    /** \name Constructors and destructors  */
    /* @{ */

    /** The default constructor */
    PerfTestVarPar(const std::string& name, PerfTestMCMC* test);

    /** The default destructor */
    ~PerfTestVarPar();

    /* @} */
    /** \name Member functions (Set)  */
    /* @{ */

    /** Set the variation parameter.
     * @param par the parameter value
     * @param name the name of the varied parameter.
     * @return an error code. */
    virtual int SetVarPar(double /*value*/, const std::string& /*name*/)
    { return 0; };

    virtual void SetProposal(bool multivariate, double dof)
    {
        fTest->SetProposal(multivariate, dof);
    }

    /* @} */
    /** \name Member functions (Get)  */
    /* @{ */

    /** Add variation parameter.
     * @param par a vector of parameter values.
     * @param name the name of the varied parameter.
     * @return an error code */
    int AddVarPar(std::vector<double> values, const std::string& name);

    /** Return the number of variation parameters. */
    int GetNVarPar()
    { return int(fVarParValues.size()); };

    /** The underlying BCModel */
    BCModel* GetModel()
    { return fTest; }

    /* @} */

    /** Run before test.
     * @return an error code. */
    int PreTest();

    /** Run the test.
     * @return an error code. */
    int RunTest();

    /** Run after test.
    * @return an error code. */
    int PostTest();

    /** Define precision settings. */
    virtual void PrecisionSettings(PerfTest::Precision);

protected:

    /** A container of varation parameters. */
    std::vector<double> fVarParValues;

    /** Name of the variation parameter. */
    std::string fVarParName;
private:

    /** The associated test. */
    PerfTestMCMC* fTest;

    /** A container of graphs. */
    std::vector<TGraphErrors*> fTargetContainer;

    /** A container of graphs. */
    std::vector<TGraphErrors*> fTestContainer;

};

#endif

