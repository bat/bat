#ifndef __BCEFFICIENCYFITTER__H
#define __BCEFFICIENCYFITTER__H

/**
 * @class BCEfficiencyFitter
 * @brief A class for fitting histograms with functions
 * @author Daniel Kollar
 * @author Kevin Kr&ouml;ninger
 * @version 1.0
 * @date 11.2008
 * @details This class allows fitting of efficiencies defined as
 * a ratio of two TH1D histograms using a TF1 function. It uses
 * binomial probabilities calculated based on the number of entries
 * in histograms. This is only applicable if the numerator is
 * a subset of the denominator.
 */

/*
 * Copyright (C) 2007-2015, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

// ---------------------------------------------------------

#include "BCFitter.h"

#include <TH1D.h>

#include <vector>

// ROOT classes
class TGraphAsymmErrors;

// ---------------------------------------------------------

class BCEfficiencyFitter : public BCFitter
{
public:

    /** \name Enumerator */
    /** @{ */

    enum DataPointType {
        kDataPointRMS = 0,              ///< Draw mean and standard deviation
        kDataPointSmallestInterval = 1, ///< Draw mean and smallest 68% interval
        kDataPointCentralInterval = 2   ///< Draw mean and central 68% interval
    };

    /** @} */

    /**Abstract class which doesn't do anything
     * but offers the right interface
     * to allow calculation the distribution of any statistic.
     * User has to create a subclass and implement the operator().
     */
    class ToyDataInterface
    {
    public:
        /**
         * operator() is called for each generated toy data set of the fast p-value calculation.
         * @param expectation the expected number of events for the parameter values
         * chosen in the call to CalculatePValueFast
         * @param toyData one toy data set */
        virtual void operator()(const std::vector<double>& expectation, const std::vector<unsigned>& toyData) = 0;

        /** pure abstract */
        virtual ~ToyDataInterface()
        { }
    };

    /** \name Constructors and destructor */
    /* @{ */

    /**
     * Constructor.
     * @param trials The histogram with the number of trials.
     * @param successes The histogram with the number of successes.
     * @param func The fit function.
     * @param name name fo the model */
    BCEfficiencyFitter(const TH1& trials, const TH1& successes, const TF1& func, const std::string& name = "efficiency_fitter_model");

    /**
     * The default destructor. */
    virtual ~BCEfficiencyFitter();

    /* @} */
    /** \name Member functions (get) */
    /* @{ */

    /**
     * @return The histogram with the number of trials */
    TH1& GetTrials()
    { return fTrials; };

    /**
     * @return The histogram with the number of successes. */
    TH1& GetSuccesses()
    { return fSuccesses; };

    /**
     * Calculates the central value and the lower and upper limits for a given probability.
     * @param n n for the binomial.
     * @param k k for the binomial.
     * @param p The central probability defining the limits.
     * @param xmin The central value.
     * @param xmin The lower limit.
     * @param xmax The upper limit.
     * @return Success of action. */
    bool GetUncertainties(int n, int k, double p, double& xexp, double& xmin, double& xmax);

    /* @} */
    /** \name Member functions (set) */
    /* @{ */

    /** Set type of point to be used to plot the efficiency data */
    void SetDataPointType(DataPointType type)
    { fDataPointType = type; }

    /* @} */
    /** \name Member functions (miscellaneous methods) */
    /* @{ */

    /**
     * The log of the prior probability. Overloaded from BCModel.
     * @param parameters A vector of doubles containing the parameter values. */
    //      virtual double LogAPrioriProbability(const std::vector<double> & parameters);

    /**
     * The log of the conditional probability. Overloaded from BCModel.
     * @param parameters A vector of doubles containing the parameter values. */
    virtual double LogLikelihood(const std::vector<double>& parameters);

    /**
     * Performs the fit.
     * @return Success of action. */
    virtual void Fit();

    /**
     * Draw the fit in the current pad. */
    virtual void DrawFit(const std::string& options = "", bool flaglegend = false);

    /**
     * Calculate the p-value using fast-MCMC. In every iteration, a new toy data set is created.
     * By providing a suitable implementation of ToyDataInterface, the user can
     * calculate the distribution of an arbitrary statistic. Each toy data set as well as the
     * expected values for the parameter values are passed on to the interface.
     * @param par A set of parameter values
     * @param pvalue The p-value for the default likelihood statistic
     * @param pvalueCorrected The p-value for the default likelihood statistic
     *                        with corrections for DoF to reduce fitting bias, in all nonpathological
     *                        cases the correction lowers the p-value.
     * @param callback requires class with operator(...) defined.
     * @param nIterations number of toy data sets generated by the Markov chain
     * @return Success of action. */
    bool CalculatePValueFast(const std::vector<double>& par, BCEfficiencyFitter::ToyDataInterface* callback, double& pvalue, double& pvalueCorrected, unsigned nIterations = 100000);

    /**
     * Calculate the p-value using fast-MCMC.
     * @param par A set of parameter values
     * @param  pvalue The pvalue
     * @param pvalueCorrected The p-value for the default likelihood statistic
     *                        with corrections for DoF to reduce fitting bias
     * @param nIterations number of pseudo experiments generated by the Markov chain
     * @return Success of action. */
    bool CalculatePValueFast(const std::vector<double>& par, double& pvalue, double& pvalueCorrected, unsigned nIterations = 100000);

    /* @} */

private:
    /**
     * The histogram containing the larger numbers. */
    TH1D fTrials;

    /**
     * The histogram containing the smaller numbers. */
    TH1D fSuccesses;

    /**
     * Temporary histogram for calculating the binomial quantiles */
    TH1D fHistogramBinomial;

    /** Type of point to plot for efficiency data */
    DataPointType fDataPointType;
};

// ---------------------------------------------------------

#endif
