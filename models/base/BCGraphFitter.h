#ifndef __BCMODELGRAPHFITTER__H
#define __BCMODELGRAPHFITTER__H

/*!
 * \class BCGraphFitter
 * \brief A class for fitting graphs with functions
 * \author Daniel Kollar
 * \author Kevin Kr&ouml;ninger
 * \version 1.0
 * \date 2008
 * \detail This class allows fitting of a TGraphErrors using
 * a TF1 function. It doeasn't take the x uncertainties into account.
 * For that look at BCGraphXFitter (not yet implemented).
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

#include <string>
#include <vector>

class TF1;
class TGraphErrors;

// ---------------------------------------------------------

class BCGraphFitter : public BCFitter
{
public:

    /** \name Constructors and destructors */
    /* @{ */

    /**
     * Constructor
     * @param name name of the model */
    BCGraphFitter(std::string name = "graph_fitter_model");

    /**
     * Constructor
     * @param graph pointer to TGraphErrors
     * @param func pointer to TF1
     * @param name name of the model */
    BCGraphFitter(TGraphErrors* graph, TF1* func, std::string name = "");

    /**
     * The default destructor. */
    virtual ~BCGraphFitter();

    /* @} */

    /** \name Member functions (get) */
    /* @{ */

    /**
     * @return pointer to TGraphErrors */
    TGraphErrors* GetGraph()
    { return fGraph; };

    /* @} */

    /** \name Member functions (set) */
    /* @{ */

    /**
     * @param graph pointer to TGraphErrors object */
    int SetGraph(TGraphErrors* graph);

    /* @} */
    /** \name Member functions (miscellaneous methods) */
    /* @{ */

    /**
     * The log of the conditional probability.
     * @param parameters vector containing the parameter values */
    double LogLikelihood(const std::vector<double>& parameters);

    /**
     * Performs the fit. The graph and the function has to be set beforehand.
     * @return Success of action. */
    bool Fit();

    /**
     * Performs the fit of the graph with the function.
     * @param graph pointer to TGraphErrors object
     * @param func pointer to TF1 object
     * @return Success of action. */
    bool Fit(TGraphErrors* graph, TF1* func);

    /**
     * Draw the fit in the current pad. */
    void DrawFit(const char* options = "", bool flaglegend = false);

    /**
     * Cumulative distribution function.
     * @param parameters
     * @param index
     * @param lower
     * @return */
    virtual double CDF(const std::vector<double>& parameters,  int index, bool lower = false);

    /**
     * Calculate chi^2, the sum of [(y-f(x))/sigma_y]^2 for all data points.
     * @param pars Parameter set to evaluate function values with.
     * @return chi^2 evalued with pars. */
    virtual double CalculateChi2(const std::vector<double>& pars);


    using BCFitter::GetPValue;

    /**
     * Calculate p value from chi^2 distribution,
     * with assumption of Gaussian distribution for all data points.
     * @param pars Parameter set to calculate p value of.
     * @param ndf Flag for choosing to include numbers of degrees of freedom.
     * @return p value if successful, negative is unsuccessful. */
    virtual double CalculatePValue(const std::vector<double>& pars, bool ndf = true);

    /* @} */

private:

    /**
     * The graph containing the data. */
    TGraphErrors* fGraph;

};

// ---------------------------------------------------------

#endif
