#ifndef __BCGOFTEST__H
#define __BCGOFTEST__H

/*!
 * \class BCGoFTest
 * \brief The class for testing model hypotheses
 * \author Daniel Kollar
 * \author Kevin Kr&ouml;ninger
 * \version 1.0
 * \date 08.2008
 * \detail This class is used for calculating the p-value of a model.
 */

/*
 * Copyright (C) 2007-2014, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

// ---------------------------------------------------------

#include "BCModel.h"

#include <string>
#include <vector>
#include <utility>

// ROOT classes
class TH1;

// BAT classes
class BCDataSet;

// ---------------------------------------------------------

class BCGoFTest : public BCModel
{
public:

    /** \name Constructors and destructors */
    /** @{ */

    /**
     * Default constructor.
     */
    BCGoFTest(std::string name);

    /**
     * Constructor with model.
     * If name is unset, "GoFTest" is appended to the model name. */
    BCGoFTest(BCModel* model, std::string name = "");

    /**
     * Default destructor. */
    ~BCGoFTest();

    /** @} */
    /** \name Member functions (get) */
    /** @{ */

    /**
     * @return distribution of log(likelihood) */
    TH1* GetHistogramLogProb()
    { return (GetNParameters() < fH1Marginalized.size()) ? fH1Marginalized[GetNParameters()] : 0; }

    /**
     * @return pointer to the tested model */
    BCModel* GetTestModel()
    { return fTestModel; };

    /**
     * @return Last p value calculated. */
    double GetPValue()
    { return fPValue; }

    /**
     * @return Last parameter set values tested. */
    std::vector<double> GetTestedParameters()
    { return fOriginalParameters; }

    /**
     * @return Log likelihood of last parameter set tested. */
    double GetTestedLogLikelihood()
    { return fOriginalLogLikelihood; }

    /** @} */
    /** \name Member functions (set) */
    /** @{ */

    /**
     * Set the model to be tested.
     * @param testmodel pointer to the model to be tested */
    void SetTestModel(BCModel* testmodel)
    { fTestModel = testmodel; };


    /** @} */
    /** \name Member functions (miscellaneous methods) */
    /** @{ */

    /**
     * Calculate the p-value.
     * @param flag_histogram A histogram is either filled or not.
     * @return p-value */
    double CalculatePValue(const std::vector<double>& parameters);

    /**
     * Updates test model's data set according to this object's parameter set,
     * then calls test model's LogLilelihood function. */
    double LogLikelihood(const std::vector<double>& parameters);

    /**
     * @return Constant value of zero. */
    double LogAPrioriProbability(const std::vector<double>& /*parameters*/)
    { return 0; };

    /**
     * Calculate the log(prob) as an observable. */
    void CalculateObservables(const std::vector<double>& pars);

    /**
     * Calls BCEngineMCMC's MCMCMetropolisPreRun(),
     * then updates boundaries of observable "log_prob"
     * for histogramming the log of the lilelihood. */
    void MCMCMetropolisPreRun();

    /** @} */

private:

    /**
     * A map of data points and data values. */
    std::vector<std::pair<int, int> > fDataMap;

    /**
     * P Value */
    double fPValue;

    /**
     * Counter for the evaluation of the p-value. */
    int fPValueBelow;
    int fPValueAbove;

    /**
     * A pointer to the model which is tested. */
    BCModel* fTestModel;

    /**
     * Original parameters to be evaluated. */
    std::vector<double> fOriginalParameters;

    /**
     * The log(likelihood) and its range. */
    double fOriginalLogLikelihood;

};

// ---------------------------------------------------------

#endif
