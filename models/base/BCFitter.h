#ifndef __BCFITTER__H
#define __BCFITTER__H

/**
 * @class BCFitter
 * @brief A base class for all fitting classes
 * @details This a general class around fitting a 1D function to data of various kinds with uncertainty propagation.
 */

/*
 * Copyright (C) 2007-2015, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

// ---------------------------------------------------------

#include "../../BAT/BCAux.h"
#include "../../BAT/BCDataSet.h"
#include "../../BAT/BCModel.h"

#include <TF1.h>
#include <TH2D.h>

#include <string>

// ---------------------------------------------------------

class BCFitter : public BCModel
{
public:

    /** \name Constructors and destructors */
    /* @{ */

    /**
     * Constructor
     * @param f pointer to TF1 to copy into fitter object
     * @param name name of the model. If empty, take the name from f. */
    BCFitter(const TF1& f, const std::string& name = "fitter_model");

    /**
     * The default destructor. */
    virtual ~BCFitter() = 0;

    /* @} */

    /** \name Member functions (get) */
    /* @{ */

    /**
     * Get fit function. */
    TF1& GetFitFunction()
    {
        return fFitFunction.at(GetCurrentChain());
    }

    /**
     * @param level Desired probability mass
     * @param nsmooth Number of times to smooth the histogram
     * @param overcoverage Flag for whether to overcover desired probability mass.
     * @return A 2D histogram of the smallest interval in Y for each bin in X containing the desired probability mass. */
    TH2* GetGraphicalErrorBandXY(double level = .68, int nsmooth = 0, bool overcoverage = true) const;

    //* @return A const reference of the internal errorband histogram
    const TH2& GetErrorBandXY()const
    {
        return fErrorBandXY;
    };

    /**
     * Returns a vector of y-values at a certain probability level.
     * @param level The level of probability
     * @return vector of y-values */
    std::vector<double> GetErrorBand(double level) const;

    TGraph* GetErrorBandGraph(double level1, double level2) const;

    TGraph* GetFitFunctionGraph(const std::vector<double>& parameters);

    TGraph* GetFitFunctionGraph()
    {
        return GetFitFunctionGraph(std::vector<double>(GetBestFitParameters().begin(), GetBestFitParameters().begin() + GetNParameters()));
    }

    TGraph* GetFitFunctionGraph(const std::vector<double>& parameters, double xmin, double xmax, int n = 1000);

    void FixDataAxis(unsigned int index, bool fixed);

    bool GetFixedDataAxis(unsigned int index) const;

    double GetPValue() const
    {
        return fPValue;
    }

    /* @} */

    /** \name Member functions (set) */
    /* @{ */

    /**
     * Sets the error band flag to continuous function */
    void SetErrorBandContinuous(bool flag);

    /**
     *Extends the lower x Edge of th errorband by -extension */
    void SetErrorBandExtensionLowEdgeX(double extension)
    {
        fErrorBandExtensionLowEdgeX = extension;
    }

    /**
     *Extends the lower x Edge of th errorband by +extension */
    void SetErrorBandExtensionUpEdgeX(double extension)
    {
        fErrorBandExtensionUpEdgeX = extension;
    }

    /**
     *Extends the lower y Edge of th errorband by -extension */
    void SetErrorBandExtensionLowEdgeY(double extension)
    {
        fErrorBandExtensionLowEdgeY = extension;
    }

    /**
     *Extends the lower y Edge of th errorband by +extension */
    void SetErrorBandExtensionUpEdgeY(double extension)
    {
        fErrorBandExtensionUpEdgeY = extension;
    }
    /**
     * Turn on or off the filling of the error band during the MCMC run.
     * @param flag set to true for turning on the filling */
    void SetFillErrorBand(bool flag = true)
    {
        fFlagFillErrorBand = flag;
    }

    /**
     * Turn off filling of the error band during the MCMC run.
     * This method is equivalent to SetFillErrorBand(false) */
    void UnsetFillErrorBand()
    {
        fFlagFillErrorBand = false;
    }

    /**
     * Sets index of the x values in function fits.
     * @param index Index of the x values */
    void SetFitFunctionIndexX(int index)
    {
        fFitFunctionIndexX = index;
    }

    /**
     * Sets index of the y values in function fits.
     * @param index Index of the y values */
    void SetFitFunctionIndexY(int index)
    {
        fFitFunctionIndexY = index;
    }

    /**
     * Sets indices of the x and y values in function fits.
     * @param indexx Index of the x values
     * @param indexy Index of the y values */
    void SetFitFunctionIndices(int indexx, int indexy)
    {
        SetFitFunctionIndexX(indexx);
        SetFitFunctionIndexY(indexy);
    }

    /**
     * Sets the flag for integration. \n
     * true: use ROOT's TF1::Integrate() \n
     * false: use linear interpolation */
    void SetFlagIntegration(bool flag)
    {
        fFlagIntegration = flag;
    };

    /**
     * Defines a fit function.
     * @param parameters A set of parameter values
     * @param x A vector of x-values
     * @return The value of the fit function at the x-values given a set of parameters */
    virtual double FitFunction(const std::vector<double>& x, const std::vector<double>& parameters);

    /**
     * Compute the integral of the fit function between xmin and xmax.
     *
     * @note Internally, the fit function needs to be integrated over the range.
     *       The integration algorithm can be toggled with SetFlagIntegration().
     */
    double Integral(const std::vector<double>& parameters, double xmin, double xmax);


    /* @} */
    /** \name Member functions (miscellaneous methods) */
    /* @{ */
    /**
     * Fit the function's parameters to the data */
    virtual void Fit() = 0;

    /**
     * Draw the fit in the current pad. */
    virtual void DrawFit(const std::string& options, bool flaglegend = false) = 0;

    /**
     * Overloaded from BCEngineMCMC */
    virtual void MCMCUserIterationInterface();

    /**
     * Overloaded from BCIntegrate. */
    virtual void MarginalizePreprocess();

    /**
     * Fill error band histogram for current iteration. This method is called from MCMCUserIterationInterface() */
    void FillErrorBand();

    /**
     * Prints a short summary of the fit results on the screen. */
    void PrintShortFitSummary();

    /* @} */

    /**
     * Create enough TF1 copies for thread safety */
    virtual void MCMCUserInitialize();

private:
    /** Fit function (as vector for thread safety) */
    std::vector<TF1> fFitFunction;

protected:
    /**
     * Copy over the most important properties of a 1D histogram.
     * This function is to emulate the TH1::Copy() method for our purposes.
     * It is not needed for root 5.34.19 and higher.
     */
    static void CopyHist(const TH1& source, TH1D& destination);

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
     * p value for goodness of fit */
    double fPValue;

    /**
     * The error band histogram */
    TH2D fErrorBandXY;

    /** Needed for uncertainty propagation */
    BCDataSet fFitterDataSet;

    /**
     * Flag for using the ROOT TH1::Integral method (true), or linear
     * interpolation (false) */
    bool fFlagIntegration;

    /** Storage for plot objects with proper clean-up */
    mutable BCAux::BCTrash<TObject> fObjectTrash;

    /**
     * extends the lower edge of x range by the given value  */
    double fErrorBandExtensionLowEdgeX;
    /**
     * extends the upper edge of x range by the given value  */
    double fErrorBandExtensionUpEdgeX;

    /**
     * extends the upper edge of y range by the given value  */
    double fErrorBandExtensionLowEdgeY;

    /**
     * extends the upper edge of y range by the given value  */
    double fErrorBandExtensionUpEdgeY;

    /** Don't allow user to accidentally set the data set,
     * as it is used internally. */
    using BCModel::SetDataSet;
    using BCModel::GetDataSet;
};

// ---------------------------------------------------------

#endif
