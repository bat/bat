#ifndef __SUMMARYTOOL__H
#define __SUMMARYTOOL__H

/*!
 * \class SummaryTool

 * This class can be used to summarize the results of an analysis. The
 * prior and posterior probabilities are compared.
 * \brief A class for summarizing the results of an analysis.
 * \author Daniel Kollar
 * \author Kevin Kr&ouml;ninger
 * \version 1.0.0
 * \date 15.02.2010
 */

/*
 * Copyright (C) 2008, 2009, 2010, Daniel Kollar and Kevin Kroeninger.
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

// ---------------------------------------------------------

#include "PriorModel.h"

#include <BAT/BCModel.h>

#include <string>

// ---------------------------------------------------------

class SummaryTool
{
 public:

	/** \name Constructors and destructors */
	/* @{ */

	/**
	 * The default constructor. */
	SummaryTool();

	/**
	 * A constructor. */
	SummaryTool(BCModel* model);

	/**
	 * The default destructor. */
	~SummaryTool();

	/* @} */
	/** \name Member functions (get) */
	/* @{ */

	/* @} */
	/** \name Member functions (set) */
	/* @{ */

	/**
	 * Set the model to be summarized.
	 * @param model The BCModel to be summarized.*/
	void SetModel(BCModel* model)
	{ fModel = model; };

	/* @} */
	/** \name Member functions (misc) */
	/* @{ */

	/**
	 * Calculate the marginalized distributions using the prior
	 * knowledge alone.
	 * @return An error flag.
	 */
	int CalculatePriorModel();

	/**
	 * Copy the summary information from the model.
	 * @return An error flag. */
	int CopySummaryData();

	/**
	 * Print a summary plot for the parameters.
	 * @return An error flag. */
	int PrintParameterPlot(const char* filename = "parameters.eps");

	/**
	 * Print a correlation plot for the parameters.
	 * @return An error flag. */
	int PrintCorrelationPlot(const char* filename = "correlation.eps");

	/**
	 * Print a comparison of the prior knowledge to the posterior
	 * knowledge for each parameter.
	 * @return An error flag. */
	int PrintKnowlegdeUpdatePlot(const char* filename = "update.eps");

	/**
	 * Print parameter summary as text.
	 * @return An error flag.*/

	int PrintParameterSummary() { return 1; };

	/**
	 * Print correlation summary as text.
	 * @return An error flag. */
	int PrintCorrelationSummary() { return 1; };

	/**
	 * Print a Latex table of the parameters.
	 * @return An error flag. */
	int PrintParameterLatex(const char* filename);

	/**
	 * Print a Latex table of the correlations.
	 * @return An error flag. */
	int PrintCorrelationLatex() { return 1; };

	/* @} */

 private:

	/**
	 * The model which results are summarized.
	 */
	BCModel* fModel;

	/**
	 * Parameter names.
	 */
	std::vector <std::string> fParName;

	/**
	 * Parameter minimum.
	 */
	std::vector <double> fParMin;

	/**
	 * Parameter maximum.
	 */
	std::vector <double> fParMax;

	/**
	 * Correlation coefficients.
	 * Length of vector equals number of parameters * number of parameters.
	 */
	std::vector <double> fCorrCoeff;

	/**
	 * Marginalized modes.\n
	 * Length of vector equals number of parameters.
	 */
	std::vector <double> fMargMode;

	/**
	 * Mean values.\n
	 * Length of vector equals number of parameters.
	 */
	std::vector <double> fMean;

	/**
	 * Global modes.\n
	 * Length of vector equals number of parameters.
	 */
	std::vector <double> fGlobalMode;

	/**
	 * Quantiles.\n
	 * The following quantiles are stored: 0.05, 0.10, 0.16, 0.5, 0.84, 0.90, 0.95.\n
	 * Length of vector equals number of parameters * number of quantiles.
	 */
	std::vector <double> fQuantiles;

	/**
	 * Smallest intervals.\n
	 * For each parameter a set of the smallest intervals is recorded.\n
	 * Structure: number of intervals n + n * (start, stop, local max, local max pos, integral)
	 * Length of vector equals number of parameters * number of quantiles.
	 */
	std::vector <double> fSmallInt;

	/**
	 * RMS values.\n
	 * Length of vector equals number of parameters.
	 */
	std::vector <double> fRMS;

	/**
	 * Sum of probabilities for quantiles
	 */
	std::vector <double> fSumProb;

	/**
	 * A model for calculating the marginalized distributions for the
	 * prior probabilities.
	 */
	PriorModel* fPriorModel;

	/**
	 * A flag: check if marginalized information is present
	 */
	bool fFlagInfoMarg;


	/**
	 * A flag: check if optimization information is present
	 */
	bool fFlagInfoOpt;

};
// ---------------------------------------------------------

#endif

