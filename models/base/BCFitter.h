#ifndef __BCFITTER__H
#define __BCFITTER__H

/*!
 * \class BCFitter
 * \brief A base class for all fitting classes
 * \author Kevin Kr&ouml;ninger
 * \version 1.0
 * \date 08.2013
 * \detail A base class for all fitting classes
 */

/*
 * Copyright (C) 2007-2014, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

// ---------------------------------------------------------

class TF1;
class TGraph;
class TH2D;

#include "../../BAT/BCModel.h"

#include <string>

// ---------------------------------------------------------

class BCFitter : public BCModel
{
public:

    /** \name Constructors and destructors */
    /* @{ */

    /**
     * Constructor
     * @param name name of the model */
    BCFitter(std::string name = "fitter_model");

    /**
     * Constructor
     * @param f pointer to TF1 to copy into fitter object
     * @param name name of the model */
    BCFitter(TF1* f, std::string name = "");

    /**
     * The default destructor. */
    virtual ~BCFitter();

    /* @} */

    /** \name Member functions (get) */
    /* @{ */

    /**
     * Get fit function
     * @param index index of chain to return fit function for (default = 0) */
    TF1* GetFitFunction(unsigned index = 0)
    { return (fFitFunction.empty()) ? (TF1*)0 : fFitFunction[index]; }

    /**
    * @return pointer to the error band */
    TGraph* GetErrorBand()
    { return fErrorBand; };

    /**
     * @return pointer to a graph for the fit function */
    TGraph* GetGraphFitFunction()
    { return fGraphFitFunction; };

    /**
     * const BCParameter * GetParameter(const char * name);
     * @return The 2D histogram of the error band. */
    TH2D* GetErrorBandXY() const
    { return fErrorBandXY; }

    /**
     * @param level Desired probability mass
     * @param nsmooth Number of times to smooth the histogram
     * @param overcoverage Flag for whether to overcover desired probability mass.
     * @return A 2D histogram of the smallest interval in Y for each bin in X containing the desired probability mass. */
    TH2D* GetGraphicalErrorBandXY(double level = .68, int nsmooth = 0, bool overcoverage = true) const;

    /**
     * Returns a vector of y-values at a certain probability level.
     * @param level The level of probability
     * @return vector of y-values */
    std::vector<double> GetErrorBand(double level) const;

    TGraph* GetErrorBandGraph(double level1, double level2) const;

    TGraph* GetFitFunctionGraph(const std::vector<double>& parameters);

    TGraph* GetFitFunctionGraph()
    { return GetFitFunctionGraph(std::vector<double>(GetGlobalMode().begin(), GetGlobalMode().begin() + GetNParameters())); }

    TGraph* GetFitFunctionGraph(const std::vector<double>& parameters, double xmin, double xmax, int n = 1000);

    void FixDataAxis(unsigned int index, bool fixed);

    bool GetFixedDataAxis(unsigned int index) const;

    double GetPValue() const
    {	return fPValue; }

    /* @} */

    /** \name Member functions (set) */
    /* @{ */

    /**
     * Set fit function
     * @param f pointer to TF1 to copy into object.
     * @return success of action. */
    bool SetFitFunction(TF1* func);

    /**
     * Turn on or off the filling of the error band during the MCMC run.
     * @param flag set to true for turning on the filling */
    void SetFillErrorBand(bool flag = true)
    { fFlagFillErrorBand = flag; }

    /**
     * Sets errorband histogram */
    void SetErrorBandHisto(TH2D* h)
    { fErrorBandXY = h; }

    /**
     * Turn off filling of the error band during the MCMC run.
     * This method is equivalent to SetFillErrorBand(false) */
    void UnsetFillErrorBand()
    { fFlagFillErrorBand = false; }

    /**
     * Sets index of the x values in function fits.
     * @param index Index of the x values */
    void SetFitFunctionIndexX(int index)
    { fFitFunctionIndexX = index; }

    /**
     * Sets index of the y values in function fits.
     * @param index Index of the y values */
    void SetFitFunctionIndexY(int index)
    { fFitFunctionIndexY = index; }

    /**
     * Sets indices of the x and y values in function fits.
     * @param indexx Index of the x values
     * @param indexy Index of the y values */
    void SetFitFunctionIndices(int indexx, int indexy)
    { SetFitFunctionIndexX(indexx); SetFitFunctionIndexY(indexy); }

    /**
     * Sets the error band flag to continuous function */
    void SetErrorBandContinuous(bool flag);

    /**
     * Defines a fit function.
     * @param parameters A set of parameter values
     * @param x A vector of x-values
     * @return The value of the fit function at the x-values given a set of parameters */
    virtual double FitFunction(const std::vector<double>& /*x*/, const std::vector<double>& /*parameters*/)
    { return 0; }

    /**
     * 1dim cumulative distribution function of the probability
     * to get the data f(x_i|param) for a single measurement, assumed to
     * be of identical functional form for all measurements
     * @param parameters The parameter values at which point to compute the cdf
     * @param index The data point index starting at 0,1...N-1
     * @param lower only needed for discrete distributions!
     * Return the CDF for the count one less than actually observed, e.g.
     * in Poisson process, if 3 actually observed, then CDF(2) is returned */
    virtual double CDF(const std::vector<double>& /*parameters*/,  int /*index*/, bool /*lower=false*/)
    {return 0.0;}

    /* @} */
    /** \name Member functions (miscellaneous methods) */
    /* @{ */

    /**
     * Read */
    int ReadErrorBandFromFile(const char* file);

    /**
     * Performs the fit.
     * @return Success of action. */
    virtual bool Fit() = 0;

    /**
     * Draw the fit in the current pad. */
    virtual void DrawFit(const char* options, bool flaglegend = false) = 0;

    /**
     * Overloaded from BCEngineMCMC */
    void MCMCIterationInterface();

    /**
     * Overloaded from BCIntegrate. */
    void MarginalizePreprocess();

    /**
     * Overloaded from BCIntegrate. */
    void MarginalizePostprocess()
    {}

    /**
     * Fill error band histogram for curreent iteration. This method is called from MCMCIterationInterface() */
    void FillErrorBand();

    /**
     * Prints a short summary of the fit results on the screen. */
    void PrintShortFitSummary();

    /* @} */

    /**
     * Sync TF1 pointers for thread safety */
    bool MCMCUserInitialize();

protected:

    /** Fit function (as vector for thread safety) */
    std::vector<TF1*> fFitFunction;

    /**
     * Pointer to a graph for displaying the fit function */
    TGraph* fGraphFitFunction;

    /**
     * Pointer to the error band (for legend) */
    TGraph* fErrorBand;

    /**
     * Flag whether or not to fill the error band */
    bool fFlagFillErrorBand;

    /**
     * The indices for function fits */
    int fFitFunctionIndexX;
    int fFitFunctionIndexY;

    /**
     * A flag for single point evaluation of the error "band" */
    bool fErrorBandContinuous;
    std::vector<double> fErrorBandX;

    /**
     * Number of X bins of the error band histogram */
    unsigned fErrorBandNbinsX;

    /**
     * Number of Y bins of the error band histogram */
    unsigned fErrorBandNbinsY;

    /**
     * The error band histogram */
    TH2D* fErrorBandXY;

    double fPValue;

};

// ---------------------------------------------------------

#endif
