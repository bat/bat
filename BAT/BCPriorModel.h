#ifndef __BCPRIORMODEL__H
#define __BCPRIORMODEL__H

/**
 * @class BCPriorModel
 * @brief Class for sampling from prior of a BCModel
 * @author Daniel Greenwald
 * @version 1.0
 * @date 09.2014
 * @details This class acts as a BCModel using the prior of another BCModel as its posterior
 * for the purpose of knowledge-update plotting.
 */

/*
 * Copyright (C) 2007-2015, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

// ---------------------------------------------------------

#include <limits>

#include "BCModel.h"

// ---------------------------------------------------------

class BCPriorModel : public BCModel
{

public:

    /** \name Constructor and destructor */
    /** @{ */

    /**
     * constructor.
     * @param model Model to be prior model of.
     * @param call_likelihood Flag to control calling of Model's likelihood. */
    BCPriorModel(BCModel& model, bool call_likelihood = false);

    /**
     * destructor. */
    virtual ~BCPriorModel() {};

    /** @} */

    /**
     * Returns a constant prior. */
    virtual double LogAPrioriProbability(const std::vector<double>& /*parameters*/)
    { return 0; }

    /**
      * Returns prior of model as posterior of PriorModel. */
    virtual double LogLikelihood(const std::vector<double>& parameters)
    { return fModel.LogAPrioriProbability(parameters); }

    /**
     * Calculates user observables according to the model. */
    void CalculateObservables(const std::vector<double>& parameters);

    /**
     * Prepare PriorModel from Model. */
    bool PreparePriorModel();

    /**
     * Set calling of likelihood in model. */
    void SetCallLikelihood(bool cl)
    { fCallLikelihood = cl; }

    /**
     * @return whether to call model's likelihood. */
    bool GetCallLikelihood() const
    { return fCallLikelihood; }

protected:
    BCModel& fModel;

    bool fCallLikelihood;

    /** \name Hide BCModel functions related to BCPriorModel. */
    /** @{ */
private:

    using BCModel::GetPriorModel;
    using BCModel::GetPrior;
    using BCModel::PrintKnowledgeUpdatePlots;

    /** @} */


};

#endif
