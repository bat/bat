#ifndef __BCMODEL__H
#define __BCMODEL__H

/**
 * @class BCModel
 * @brief The base class for all user-defined models.
 * @author Daniel Kollar
 * @author Kevin Kr&ouml;ninger
 * @author Daniel Greenwald
 * @version 1.0
 * @date 08.2008
 * @details This class represents a model. It contains a container of prior distributions and the likelihood. The methods that implement the prior and the likelihood
 * have to be overloaded by the user in the user defined model class
 * derived from this class.
 */

/*
 * Copyright (C) 2007-2015, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

// ---------------------------------------------------------

#include "BCH1D.h"
#include "BCH2D.h"
#include "BCIntegrate.h"
#include "BCDataSet.h"

#include <string>

//BAT classes
class BCPriorModel;

// ---------------------------------------------------------

class BCModel : public BCIntegrate
{

public:

    /** \name Constructors and destructors */
    /** @{ */

    /**
     * Default constructor.
     * @param name The name of the model */
    BCModel(const std::string& name = "model");

    /**
     * Copy constructor. */
    BCModel(const BCModel& bcmodel);

    /**
     * Read in MCMC constructor.
     * @param filename Path of file holding model.
     * @param name Name of model (file should contain TTree's [name]_mcmc and [name]_parameters.\n
     * if empty string is given, properly matching TTrees are searched for in the file.
     * @param loadObservables Flag for whether to load observables for file (true; default) or to let user respecify observables.*/
    BCModel(const std::string& filename, const std::string& name, bool loadObservables = true);

    /**
     * Destructor. */
    virtual ~BCModel();

    /** @} */
    /** \name swap*/
    /** @{ */

    /** swap */
    friend void swap(BCModel& A, BCModel& B);

    /** @} */
    /** \name Member functions (get) */
    /** @{ */

    /**
     * @return The data set. */
    BCDataSet* GetDataSet()
    { return fDataSet; }

    /**
     * @return The number of data points in the current data set. */
    unsigned GetNDataPoints() const
    { return (fDataSet) ? fDataSet->GetNDataPoints() : 0; }

    /**
     * @return The number of degrees of freedom:
     * (number of data points) - (number of free parameters) */
    int GetNDoF() const
    { return GetNDataPoints() - fParameters.GetNFreeParameters(); }

    /**
     * @return BCPriorModel. */
    virtual BCPriorModel* GetPriorModel(bool prepare = true, bool call_likelihood = false);

    /**
     * Get prior of a variable as a BCH1D.
     * If not using factorized priors, first call GetPriorModel()->Marginalize() first.
     * @param index Index of variable to get prior for.
     * @return BCH1D of prior */
    virtual BCH1D GetPrior(unsigned index);

    /**
     * Get prior of a pair of variables as a BCH2D.
     * If not using factorized priors, first call GetPriorModel()->Marginalize() first.
     * @param index1 Index of variable to get prior for (abcissa).
     * @param index2 Index of variable to get prior for (ordinate).
     * @return BCH2D of prior */
    virtual BCH2D GetPrior(unsigned index1, unsigned index2);

    /**
     * @return BCH1D object for controlling drawing options of priors in knowledge update plots. */
    BCH1D& GetBCH1DPriorDrawingOptions()
    { return fBCH1DPriorDrawingOptions; }

    /**
     * @return BCH2D object for controlling drawing options of priors in knowledge update plots. */
    BCH2D& GetBCH2DPriorDrawingOptions()
    { return fBCH2DPriorDrawingOptions; }

    /**
     * @return BCH1D object for controlling drawing options of posteriors in knowledge update plots. */
    BCH1D& GetBCH1DPosteriorDrawingOptions()
    { return fBCH1DPosteriorDrawingOptions; }

    /**
     * @return BCH2D object for controlling drawing options of posteriors in knowledge update plots. */
    BCH2D& GetBCH2DPosteriorDrawingOptions()
    { return fBCH2DPosteriorDrawingOptions; }

    /**
     * @return flag for drawing prior first (true) or
     * posterior first (false) in knowledge update plots. */
    bool GetDrawPriorFirst() const
    { return fDrawPriorFirst; }

    /** @} */

    /** \name Member functions (set) */
    /** @{ */

    /**
     * Sets the data set. Model does not own data set!
     * @param dataset A data set */
    void SetDataSet(BCDataSet* dataset)
    { fDataSet = dataset; }

    /**
     * Set default drawing options for knowledge update plots. */
    void SetKnowledgeUpdateDrawingStyle(BCAux::BCKnowledgeUpdateDrawingStyle style = BCAux::kKnowledgeUpdateDefaultStyle);

    /**
     * Set drawing of prior first (true) or posterior first (false)
     * for knowledge update plots. */
    void SetDrawPriorFirst(bool b = true)
    { fDrawPriorFirst = b; }

    /** @} */

    /** \name Member functions (miscellaneous methods) */
    /** @{ */

    /**
     * Returns the prior probability.
     * @param parameters A set of parameter values
     * @return The prior probability p(parameters)
     * @see GetPrior(std::vector<double> parameters) */
    virtual double APrioriProbability(const std::vector<double>& parameters)
    { return exp(this->LogAPrioriProbability(parameters)); }

    /**
     * Returns natural logarithm of the prior probability.
     * Assumes prior factorizes into independent parts set by user for each parameter.
     * Method can be overloaded by the user (especially for nonfactorizable priors)
     * @param parameters A set of parameter values
     * @return The prior probability p(parameters)
     * @see GetPrior(std::vector<double> parameters) */
    virtual double LogAPrioriProbability(const std::vector<double>& parameters)
    { fFactorizedPrior = true; return fParameters.GetLogPrior(parameters); }

    /**
     * Returns the likelihood
     * @param params A set of parameter values
     * @return The likelihood */
    virtual double Likelihood(const std::vector<double>& params)
    { return exp(LogLikelihood(params)) ; }

    /**
     * Calculates natural logarithm of the likelihood.
     * Method needs to be overloaded by the user.
     * @param params A set of parameter values
     * @return Natural logarithm of the likelihood */
    virtual double LogLikelihood(const std::vector<double>& params) = 0;

    /**
     * Returns the likelihood times prior probability given a set of parameter values
     * @param params A set of parameter values
     * @return The likelihood times prior probability */
    virtual double ProbabilityNN(const std::vector<double>& params)
    { return exp(LogProbabilityNN(params)); }

    /**
     * Returns the natural logarithm of likelihood times prior probability given
     * a set of parameter values
     * @param parameters A set of parameter values
     * @return The likelihood times prior probability */
    virtual double LogProbabilityNN(const std::vector<double>& parameters);

    /**
     * Returns the a posteriori probability given a set of parameter values
     * @param parameters A set of parameter values
     * @return The posterior probability */
    virtual double Probability(const std::vector<double>& parameters)
    { return exp(LogProbability(parameter)); }

    /**
     * Returns natural logarithm of the  a posteriori probability given a set of parameter values
     * @param parameters A set of parameter values
     * @return log(posterior probability) */
    virtual double LogProbability(const std::vector<double>& parameters);

    /**
     * Sampling function used for importance sampling.
     * Method needs to be overloaded by the user.
     * @param parameters A set of parameter values
     * @return The probability density at the parameter values */
    virtual double SamplingFunction(const std::vector<double>& parameters);

    /**
     * Overloaded function to evaluate integral. */
    virtual double Eval(const std::vector<double>& parameters)
    { return exp(LogEval(parameters)); }

    /**
     * Overloaded function to evaluate integral. */
    virtual double LogEval(const std::vector<double>& parameters)
    { return LogProbabilityNN(parameters); }

    /**
     * Initialize the trees containing the Markov chains and parameter info. */
    virtual void InitializeMarkovChainTree(bool replacetree = false, bool replacefile = false);

    /**
     * Calculates the matrix element of the Hessian matrix
     * @param index1 Index of parameter for the first derivative
     * @param index2 Index of parameter for the second derivative
     * @param point Vector of parameters for point at which to evaluate.
     * @return The matrix element of the Hessian matrix */
    double HessianMatrixElement(unsigned index1, unsigned index2, const std::vector<double>& point);

    /**
     * Prints a short summary of the fit results on the screen. */
    void PrintShortFitSummary();

    /**
     * Prints matrix elements of the Hessian matrix
     * @param parameters The parameter values at which point to evaluate the matrix */
    void PrintHessianMatrix(std::vector<double> parameters);

    /**
     * Print a comparison of the prior knowledge to the posterior
     * knowledge for each parameter.
     * @param filename name of file to print to
     * @param hdiv number of horizontal divisions on page
     * @param vdiv number of vertical divisions on page
     * @param call_likelihood whether to call model likelihood function to prepare for calculating observables
     * @return number of plots printed. */
    virtual unsigned PrintKnowledgeUpdatePlots(const std::string& filename, unsigned hdiv = 1, unsigned vdiv = 1, bool call_likelihood = false);


    /* /\** */
    /*  * Print priors to file */
    /*  * @param filename name of file to print to */
    /*  * @param hdiv number of horizontal divisions on page */
    /*  * @param vdiv number of vertical divisions on page */
    /*  * @param call_likelihood whether to call model likelihood function to prepare for calculating observables */
    /*  * @return number of plots printed. *\/ */
    /* virtual unsigned PrintPriors(const std::string& filename, unsigned hdiv = 1, unsigned vdiv = 1, bool call_likelihood = false); */

    /** @} */

    // declare BCPriorModel a friend
    friend class BCPriorModel;

protected:

    /**
     * A data set */
    BCDataSet* fDataSet;

    /**
     * BCPriorModel object for drawing of knowledge update, and saving of samples according to prior.*/
    BCPriorModel* fPriorModel;

    /**
     * knowledge update plot 1D prior options. */
    BCH1D fBCH1DPriorDrawingOptions;

    /**
     * knowledge update plot 2D prior options. */
    BCH2D fBCH2DPriorDrawingOptions;

    /**
     * knowledge update plot 1D posterior options. */
    BCH1D fBCH1DPosteriorDrawingOptions;

    /**
     * knowledge update plot 2D posterior options. */
    BCH2D fBCH2DPosteriorDrawingOptions;

    /**
     * flag for ordering of drawing of prior and posterior in knowledge update plots. */
    bool fDrawPriorFirst;

    /**
     * flag for whether factorized prior has been used. */
    bool fFactorizedPrior;

};

// ---------------------------------------------------------

#endif
