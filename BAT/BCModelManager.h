#ifndef __BCMODELMANAGER__H
#define __BCMODELMANAGER__H

/**
 * @class BCModelManager
 * @brief A class representing a set of BCModels.
 * @author Daniel Kollar
 * @author Kevin Kr&ouml;ninger
 * @version 1.0
 * @date 08.2008
 * @details This class represents a manager for BCModels. It handles
 * common data sets and performs operations on BCModels
 * simultaneously. Model comparsion in terms of a posteriori
 * probabilities is only possible with this class.
 */

/*
 * Copyright (C) 2007-2015, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

// ---------------------------------------------------------

#include "BCDataSet.h"
#include "BCModel.h"

#include <string>

// ---------------------------------------------------------

class BCModelManager
{
public:

    /** \name Constructors and destructor */
    /** @{ */

    /**
     * The default constructor. */
    BCModelManager();

    /**
     * The default copy constructor. */
    BCModelManager(const BCModelManager& other);

    /**
     * The default destructor. */
    virtual ~BCModelManager();

    /** @} */

    /** \name operators and swap*/
    /** @{ */

    /**
     * The defaut assignment operator */
    BCModelManager& operator=(BCModelManager modelmanager);

    /** swap */
    friend void swap(BCModelManager& A, BCModelManager& B);

    /** @} */

    /** \name Member functions (get) */
    /** @{ */

    /**
     * @return The number of models. */
    unsigned int GetNModels()
    { return fModels.size(); }

    /**
     * Returns the BCModel at a certain index of this BCModelManager.
     * @param index The index of the model in the BCModelManager.
     * @return The BCModel at the index. */
    BCModel* GetModel(unsigned index)
    { return (index < fModels.size()) ? fModels[index] : 0; }

    /**
     * Returns the common data set.
     * @return The data set. */
    BCDataSet* GetDataSet()
    { return fDataSet; };

    /** @} */

    /** \name Member functions (set) */
    /** @{ */

    /**
     * Sets the data set common to all BCModels in this
     * BCModelManager. Note: Data set is _not_ owned by the manager,
     * nor the models.
     * @param dataset A data set */
    void SetDataSet(BCDataSet* dataset);

    /**
     * Set the precision for the MCMC run. */
    void SetPrecision(BCEngineMCMC::Precision precision);

    /**
     * Sets the maximum number of iterations for the Monte Carlo
     * integration for all BCModels in this BCModelManger. */
    void SetNIterationsMax(int niterations);

    /**
     * Sets the minimum number of iterations for the Monte Carlo
     * integration for all BCModels in this BCModelManager.
     * @param niterations */
    void SetNIterationsMin(int niterations);

    /**
     * @param niterations interval for checking precision in integration routines */
    void SetNIterationsPrecisionCheck(int niterations);

    /**
     * @param method The marginalization method */
    void SetMarginalizationMethod(BCIntegrate::BCMarginalizationMethod method);

    /**
     * @param method The integration method */
    void SetIntegrationMethod(BCIntegrate::BCIntegrationMethod method);

    /**
     * @param method The mode finding method */
    void SetOptimizationMethod(BCIntegrate::BCOptimizationMethod method);

    /**
     * @param relprecision The relative precision envisioned for Monte
     * Carlo integration */
    void SetRelativePrecision(double relprecision);

    /**
     * Set absolute precision of the numerical integation */
    void SetAbsolutePrecision(double absprecision);

    /**
     * @param n Number of bins per dimension for the marginalized
     * distributions.  Default is 100. Minimum number allowed is 2. */
    void SetNbins(unsigned int n);

    /**
     * Sets the number of Markov chains */
    void SetNChains(unsigned int n);

    /** @} */

    /** \name Member functions (miscellaneous methods) */
    /** @{ */

    /**
     * Adds a model to the container
     * @param model The model
     * @param probability The a priori probability
     * @see AddModel(BCModel * model)
     * @see SetModelPrior(BCModel * model, double probability) */
    void AddModel(BCModel* model, double prior_probability = 0.);

    /**
     * Calculates the normalization of the likelihood for each model in
     * the container. */
    void Integrate();

    /**
     * Calculate Bayes factor for two models.
     * @param imodel1 index of model 1 (numerator)
     * @param imodel2 index of model 2 (denominator)
     * @return Bayes factor or -1. on error */
    double BayesFactor(const unsigned int imodel1, const unsigned int imodel2) const;

    /**
     * Does the mode finding */
    void FindMode();

    /**
     * Marginalize all probabilities wrt. single parameters and all
     * combinations of two parameters for all models. */
    void MarginalizeAll();

    /**
     * Turn on/off writing of Markov chains to root files for all models.
     * If setting true, you must first set filename with function with filename arguments.
     * @param flag Flag for writing Markov chain (run and prerun) to ROOT file. */
    void WriteMarkovChain(bool flag);

    /**
      * Turn on/off writing of Markov chain to root file during run for all models.
      * If setting either true, you must first set filename with function with filename arguments.
      * @param flag Flag for writing run Markov chain to ROOT file (true) or not (false). */
    void WriteMarkovChainRun(bool flag);

    /**
     * Turn on/off writing of Markov chain to root file during prerun for all models.
     * If setting either true, you must first set filename with function with filename arguments.
     * @param flag Flag for writing prerun Markov chain to ROOT file (true) or not (false). */
    void WriteMarkovChainPreRun(bool flag);

    /**
     * Turn on writing of Markov chains to root files for all models.
     * @param prefix MCMC's written to files named "[prefix][Model safe name].root"
     * @param option file-open options (TFile), must be "NEW", "CREATE", "RECREATE", or "UPDATE" (i.e. writeable).
     * @param flag_run Flag for writing run Markov chain to ROOT file (true) or not (false).
     * @param flag prerun Flag for writing prerun Markov chain to ROOT file (true) or not (false). */
    void WriteMarkovChain(const std::string& filename, const std::string& option, bool flag_run = true, bool flag_prerun = true);

    /**
     * Prints a summary of the model comparison to the log. */
    void PrintModelComparisonSummary() const;

    /**
     * Prints a summary to the logs. */
    void PrintSummary() const;

    /** @} */

private:

    /**
     * Vector of pointers to all models. */
    std::vector<BCModel*> fModels;

    /**
     * Vector of a priori probabilities. */
    std::vector<double> fAPrioriProbability;

    /**
     * Vector of a posteriori probabilities. */
    std::vector<double> fAPosterioriProbability;

    /**
     * The data set common to all models. */
    BCDataSet* fDataSet;

};

// ---------------------------------------------------------

#endif
