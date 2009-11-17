#ifndef __STACKMODELMANAGER__H
#define __STACKMODELMANAGER__H

/*!
 * \class StackModelManager
 * This class can be used for managing several StackModel objects. A
 * variety of methods to compare the different models is available.
 * \brief A class for managing several StackModel objects.
 * \author Kevin Kr&ouml;ninger
 * \version 1.0
 * \date 10.2009
 */

/*
 * Copyright (C) 2008, Daniel Kollar and Kevin Kroeninger.
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

// ---------------------------------------------------------

#include <StackModel.h>

#include <TH1D.h>

// ---------------------------------------------------------

class StackModelManager
{
 public:

	/**
	 * The default constructor.
	 */
	StackModelManager();

	/**
	 * The default destructor.
	 */
	~StackModelManager();

	/**
	 * Set the histogram containing the data.
	 */
	int SetDataHistogram(TH1D hist);

	/**
	 * Return a stack model.
	 * @param index The index of the model.
	 */
	StackModel * GetStackModel(int index)
		{ return fStackModelContainer.at(index); };

	/**
	 * Add a stack model to the container.
	 */
	void AddStackModel(StackModel * sm);

	/**
	 * Return the number of models in the container.
	 */
	int GetNModels()
	{ return int(fStackModelContainer.size()); };

	/**
	 * Return the chi2 for a model.
	 */
	double GetChi2(int index)
	{ return fChi2.at(index); };

	/**
	 * Return the chi2-prob. for a model.
	 */
	double GetChi2Prob(int index)
	{ return fChi2Prob.at(index); };

	/**
	 * Return the Kolomogorov-Smirnov probability for a model.
	 */
	double GetKSProb(int index)
	{ return fKSProb.at(index); };

	/**
	 * Return the p-value for a model.
	 */
	double GetPValue(int index)
	{ return fPValue.at(index); };

	/**
	 * Return the maximum Likelihood for a model.
	 */
	double GetMaxLike(int index)
	{ return fMaxLike.at(index); };

	/**
	 * Return the posterior probabiltiy for a model for a model.
	 */
	double GetPosterior(int index)
	{ return fPosterior.at(index); };

	/**
	 * Return the normalization for a model.
	 */
	double GetNorm(int index)
	{ return fNorm.at(index); };

	/**
	 * Return the K-factor for two models
	 */
	double GetKFactor(int index1, int index2)
	{ return (fKFactor.at(index2)).at(index2); };

	/**
	 * Return the number of parameters of a model.
	 */
	int GetNParameters(int index)
	{ return fNParameters.at(index); };

	/**
	 * Return the number of degrees of freedom of a model.
	 */
	int GetNDF(int index)
	{ return fNDF.at(index); };

	/**
	 * Performs a series of analyses for the models in the container.
	 */
	int PerformAnalysis();

	/**
	 * Prints the results to screen.
	 */
	void PrintResults();

 protected:
	/**
	 * A container for the stack models.
	 */
	std::vector <StackModel *> fStackModelContainer;

	/**
	 * A vector containing the chi2 values.
	 */
	std::vector <double> fChi2;

	/**
	 * A vector containing the chi2-prob. values.
	 */
	std::vector <double> fChi2Prob;

	/**
	 * A vector containing the Kolmogorov-Smirnov probability values.
	 */
	std::vector <double> fKSProb;

	/**
	 * A vector containing the p-values.
	 */
	std::vector <double> fPValue;

	/**
	 * A vector containing the maximum Likelihood values.
	 */
	std::vector <double> fMaxLike;

	/**
	 * A vector containing the posterior probability values.
	 */
	std::vector <double> fPosterior;

	/**
	 * A vector containing the normalization factors.
	 */
	std::vector <double> fNorm;

	/**
	 * A vector containing a vector of K-factors.
	 */
	std::vector <std::vector<double> > fKFactor;

	/**
	 * A vector containing the number of parameters values.
	 */
	std::vector <int> fNParameters;

	/**
	 * A vector containing the number of degrees of freedom.
	 */
	std::vector <int> fNDF;

};
// ---------------------------------------------------------

#endif

