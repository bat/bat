#ifndef __BCMODELMANAGER__H
#define __BCMODELMANAGER__H

/*!
 * \class BCModelManager
 * \brief A class representing a set of BCModels.
 * \author Daniel Kollar
 * \author Kevin Kr&ouml;ninger
 * \version 1.0
 * \date 08.2008
 * \detail This class represents a manager for BCModels. It handles
 * common data sets and performs operations on BCModels
 * simultaneously. Model comparsion in terms of a posteriori
 * probabilities is only possible with this class.
 */

/*
 * Copyright (C) 2007-2013, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

// ---------------------------------------------------------

#include "BCModel.h"
#include "BCDataSet.h"

#include <string>

// ---------------------------------------------------------

class BCModelManager {
public:

	/** \name Constructors and destructors */
	/** @{ */

	/**
	 * The default constructor. */
	BCModelManager();

	/**
	 * The default copy constructor. */
	BCModelManager(const BCModelManager & other);

	/**
	 * The default destructor. */
	virtual ~BCModelManager();

	/** @} */

	/** \name Assignment operators */
	/** @{ */

	/**
	 * The defaut assignment operator */
	BCModelManager & operator = (const BCModelManager & modelmanager);

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
	BCModel * GetModel(unsigned index)
	{ return (index<fModels.size()) ? fModels[index] : 0; }

	/**
	 * Returns the common data set.
	 * @return The data set. */
	BCDataSet * GetDataSet()
	{ return fDataSet; };

	/** @} */

	/** \name Member functions (set) */
	/** @{ */

	/**
	 * Sets the data set common to all BCModels in this
	 * BCModelManager.
	 * @param dataset A data set */
	void SetDataSet(BCDataSet * dataset);

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
	void AddModel(BCModel * model, double probability=0.);

	/**
	 * Calculates the normalization of the likelihood for each model in
	 * the container. */
	void Integrate();

	/**
	 * Calculate Bayes factor for two models.
	 * @param imodel1 index of model 1 (numerator)
	 * @param imodel2 index of model 2 (denominator)
	 * @return Bayes factor or -1. on error */
	double BayesFactor(const unsigned int imodel1, const unsigned int imodel2);

	/**
	 * Does the mode finding */
	void FindMode();

	/**
	 * Marginalize all probabilities wrt. single parameters and all
	 * combinations of two parameters for all models. */
	void MarginalizeAll();

	/**
	 * Flag for writing Markov chain to file */
	void WriteMarkovChain(bool flag);

	/**
	 * Turn on writing of Markov chain to file for all models.
	 * @param prefix Prefix for output files (model safe name is appended along with .root)
	 * @param option file-opn options (TFile), must be "NEW", "CREATE", "RECREATE", or "UPDATE". */
	void WriteMarkovChain(std::string prefix, std::string option);

	/**
	 * Prints a summary of model comparison into a file.
	 * If filename is omitted the summary will be printed onto the screen
	 * @param filename name of the file to write into. */
	void PrintModelComparisonSummary(const char * filename=0);

	/**
	 * Prints a summary into a file. If filename is omitted the summary
	 * will be printed onto the screen.
	 * This method is obsolete. Use PrintResults() instead.
	 * @param filename name of the file to write into. */
	void PrintSummary(const char * filename=0);

	/**
	 * Prints summaries of all files */
	void PrintResults();

	/* /\** */
	/*  * Calculates the p-value for all models. *\/ */
	/* void CalculatePValue(bool flag_histogram=false); */

	/** @} */

private:

	/**
	 * Vector of pointers to all models. */
	std::vector<BCModel*> fModels;

	/**
	 * The data set common to all models. */
	BCDataSet * fDataSet;

};

// ---------------------------------------------------------

#endif
